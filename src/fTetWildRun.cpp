#include "fTetWildRun.h"

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/MshLoader.h>

#include <Eigen/Dense>
#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh_reorder.h>

using namespace floatTetWild;
using namespace Eigen;





/*
namespace floatTetWild {
uint32_t HashProbe(Mesh &mesh) {

    //from http://isthe.com/chongo/src/fnv/hash_32.c
    auto fnv_32_buf = [](const void *buf, size_t len, uint32_t hval = 0){
        const unsigned char *bp = (const unsigned char *)buf;
        const unsigned char *be = bp + len;

        while (bp < be) {
	        hval *= 0x01000193;
	        hval ^= 1 + (uint32_t)*bp++;
        }

        return hval;
    };

    //Hash a single bool. Separate because false is guaranteeed zero and true is only guaranteed nonzero
    auto fnv_32_bool = [&fnv_32_buf](bool b, uint32_t hval = 0){
        unsigned char buf = 0;
        if(b) buf = 1;

        return fnv_32_buf(&buf, 1, hval);
    };


    uint32_t hash = 0;

    //mesh
    for(const auto &e : mesh.tet_vertices){
        hash = fnv_32_buf(e.pos.data(), 3*sizeof(Scalar), hash);
        hash = fnv_32_buf(e.conn_tets.data(), e.conn_tets.size()*sizeof(int), hash);

        hash = fnv_32_bool(e.is_on_surface , hash);
        hash = fnv_32_bool(e.is_on_boundary, hash);
        hash = fnv_32_bool(e.is_on_cut     , hash);
        hash = fnv_32_bool(e.is_on_bbox    , hash);
        hash = fnv_32_bool(e.is_outside    , hash);
        hash = fnv_32_bool(e.is_removed    , hash);
        hash = fnv_32_bool(e.is_freezed    , hash);

        hash = fnv_32_buf(&e.on_boundary_e_id, sizeof(int), hash);
        hash = fnv_32_buf(&e.sizing_scalar, sizeof(Scalar), hash);
        hash = fnv_32_buf(&e.scalar, sizeof(Scalar), hash);
    }

    for(const auto &e : mesh.tets){
        hash = fnv_32_buf(e.indices.data(), 4*sizeof(int), hash);

        hash = fnv_32_buf(e.is_surface_fs.data(), 4*sizeof(char), hash);
        hash = fnv_32_buf(e.is_bbox_fs   .data(),    4*sizeof(char), hash);
        hash = fnv_32_buf(e.opp_t_ids    .data(),     4*sizeof(int), hash);
        hash = fnv_32_buf(e.surface_tags .data(),  4*sizeof(char), hash);

        hash = fnv_32_buf(&e.quality, sizeof(Scalar), hash);
        hash = fnv_32_buf(&e.scalar, sizeof(Scalar), hash);

        hash = fnv_32_bool(e.is_removed, hash);
        hash = fnv_32_bool(e.is_outside, hash);
    }

    hash = fnv_32_buf(&mesh.t_empty_start, sizeof(int), hash);
    hash = fnv_32_buf(&mesh.v_empty_start, sizeof(int), hash);
    hash = fnv_32_bool(mesh.is_limit_length      , hash);
    hash = fnv_32_bool(mesh.is_closed            , hash);
    hash = fnv_32_bool(mesh.is_input_all_inserted, hash);
    hash = fnv_32_bool(mesh.is_coarsening        , hash);

    return hash;
}
}*/







std::pair<Eigen::MatrixXf, Eigen::Matrix4Xi> fTetWildRun(const Eigen::MatrixXf &Vertices, const Eigen::Matrix3Xi &Faces, floatTetWild::Scalar ideal_edge_length_rel, std::string OutputName, bool Quiet, unsigned int max_threads){
#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    //Note: this is necessary and results in printing, at least in debug mode :(
    //It's inconsistent.
    GEO::initialize();
    bool skip_simplify = false;

    Mesh        mesh;
    Parameters& params = mesh.params;
Quiet=true;
    if(Quiet){
        params.log_level = 6;
        params.is_quiet = true;
        spdlog::set_level(spdlog::level::off);
    }

    params.ideal_edge_length_rel = ideal_edge_length_rel;

#ifdef FLOAT_TETWILD_USE_TBB
    const size_t MB          = 1024 * 1024;
    const size_t stack_size  = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(max_threads, num_threads);
    params.num_threads       = num_threads;

    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif


    /// set sizing field
    Eigen::VectorXd V_in;
    Eigen::VectorXi T_in;
    Eigen::VectorXd values;

    if (V_in.rows() != 0 && T_in.rows() != 0 && values.rows() != 0) {
        params.apply_sizing_field = true;

        params.V_sizing_field = V_in;
        params.T_sizing_field = T_in;
        params.values_sizing_field = values;
    }

    /// set input tage
    std::vector<Vector3>  input_vertices(Vertices.cols());
    std::vector<Vector3i> input_faces   (   Faces.cols());
    std::vector<int>      input_tags;

    for(size_t i = 0; i != Vertices.cols(); i++) input_vertices[i] = Vertices.block(0, i, 3, 1).cast<floatTetWild::Scalar>();
    for(size_t i = 0; i !=    Faces.cols(); i++)    input_faces[i] = Faces   .block(0, i, 3, 1);


    /// set envelope
    GEO::Mesh                sf_mesh;

    //I hate this: create a GEO::Mesh alongside.
    sf_mesh.vertices.create_vertices(Vertices.cols());
    sf_mesh.facets.create_facets(Faces.cols(), 3);

    for(size_t i = 0; i != Vertices.cols(); i++){
        auto coords = Vertices.col(i).data();

        //Content of set_mesh_point copied from GEO mesh_io
        if(sf_mesh.vertices.single_precision()) {
            float* p = sf_mesh.vertices.single_precision_point_ptr(i);
            for(size_t c=0; c<3; ++c) {
                p[c] = coords[c];
            }
        } else {
            double* p = sf_mesh.vertices.point_ptr(i);
            for(size_t c=0; c<3; ++c) {
                p[c] = double(coords[c]);
            }
        }
    }

    for(size_t i = 0; i != Faces.cols(); i++){
        auto f = Faces.col(i).data();
        sf_mesh.facets.set_vertex(i, 0, f[0]);
        sf_mesh.facets.set_vertex(i, 1, f[1]);
        sf_mesh.facets.set_vertex(i, 2, f[2]);
    }

    GEO::mesh_reorder(sf_mesh, GEO::MESH_ORDER_MORTON);


    if (input_tags.size() != input_faces.size()) {
        input_tags.resize(input_faces.size());
        std::fill(input_tags.begin(), input_tags.end(), 0);
    }


    AABBWrapper tree(sf_mesh);

    if (!params.init(tree.get_sf_diag()))
        return std::make_pair(MatrixXf(), Matrix4Xi());

    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);

    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
//cout << HashProbe(mesh) << "\n";
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
//Ben S note: does this lib contain a Delaunay tetrahedralize function? If so and if it's any good, consider just running it on a Poisson disk sampled point cloud.
//cout << HashProbe(mesh) << "\n"; exit(0);
//cout << mesh.get_avg_energy() << "\n";
//cout << mesh.get_max_energy() << "\n";
//cout << mesh.get_t_num() << "\n";
//cout << mesh.get_v_num() << "\n\n";

    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);

//cout << mesh.get_avg_energy() << "\n";
//cout << mesh.get_max_energy() << "\n";
//cout << mesh.get_t_num() << "\n";
//cout << mesh.get_v_num() << "\n\n";

    optimization(
      input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});

//cout << mesh.get_avg_energy() << "\n";
//cout << mesh.get_max_energy() << "\n";
//cout << mesh.get_t_num() << "\n";
//cout << mesh.get_v_num() << "\n\n";

    correct_tracked_surface_orientation(mesh, tree);
//cout << mesh.get_avg_energy() << "\n";
//cout << mesh.get_max_energy() << "\n";
//cout << mesh.get_t_num() << "\n";
//cout << mesh.get_v_num() << "\n\n";

    if (params.smooth_open_boundary) {
        smooth_open_boundary(mesh, tree);
        for (auto& t : mesh.tets) {
            if (t.is_outside)
                t.is_removed = true;
        }
    }
    else {
        if (!params.disable_filtering) {
            if (params.use_floodfill) {
                filter_outside_floodfill(mesh);
            }
            else if (params.use_input_for_wn) {
                filter_outside(mesh, input_vertices, input_faces);
            }
            else
                filter_outside(mesh);
        }
    }


    if(!OutputName.empty()){
        //Write tetrahedral mesh.
        MeshIO::write_mesh(OutputName + ".msh", mesh, false);

        //This computes/writes boundary surface of tetrahedral mesh
        Eigen::MatrixXd V_sf;
        Eigen::MatrixXi F_sf;
        if(mesh.params.manifold_surface)
            manifold_surface(mesh, V_sf, F_sf);
        else
            get_surface(mesh, V_sf, F_sf);

        igl::write_triangle_mesh(OutputName + "_boundary.obj", V_sf, F_sf);
    }



    //Convert to more primitive data.
    MatrixXf TV;
    Matrix4Xi TF;

    TV.resize(3, mesh.tet_vertices.size());
    TF.resize(NoChange, mesh.tets.size());

    for(size_t i = 0; i != mesh.tet_vertices.size(); i++){
        const auto &v = mesh.tet_vertices[i];
        auto p = TV.col(i).data();

        p[0] = (float)v[0];
        p[1] = (float)v[1];
        p[2] = (float)v[2];
    }

    for(size_t i = 0; i != mesh.tets.size(); i++){
        const auto &t = mesh.tets[i];
        auto p = TF.col(i).data();

        p[0] = t[0];
        p[1] = t[1];
        p[2] = t[2];
        p[3] = t[3];
    }


    return std::make_pair(TV, TF);
}


