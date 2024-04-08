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


/* As of this writing, fTetWildRun is deterministic (same output with same input) if max_threads = 1.
Todo:
    See if you can make it deterministic with more than one thread.
    Enable float (vs double). Others have done it.

Other forks of interest:

	https://github.com/Mathias-Fuchs/fTetWild
	https://github.com/seeeagull/fTetWild/tree/master

There's a build option - float vs double. Any disadvantage to float, other than it doesn't freaking work?

*/


//Only after writing this I realized FloatTetwild.cpp is basically a copy of this.

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

    if(Quiet){
        params.log_level = 6;
        params.is_quiet = true;
        logger().set_level(spdlog::level::off);
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
    std::vector<floatTetWild::Vector3>  input_vertices(Vertices.cols());
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

    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);


    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);

    optimization(
      input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});

    correct_tracked_surface_orientation(mesh, tree);

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
        std::vector<Scalar> colors;

        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }

        MeshIO::write_mesh(OutputName + ".msh", mesh, false, colors);


        //This computes/writes boundary surface of tetrahedral mesh
        Eigen::MatrixXd V_sf;
        Eigen::MatrixXi F_sf;
        if(mesh.params.manifold_surface)
            manifold_surface(mesh, V_sf, F_sf);
        else
            get_surface(mesh, V_sf, F_sf);

        igl::write_triangle_mesh(OutputName + "_boundary.obj", V_sf, F_sf);
    }


    //Convert to more primitive data. It's a bit confounded by the presense of vertices and tetrahedra that are to be removed,
    //see write_mesh_aux in MeshIO.cpp for how that's originally done. Start by counting and allocating.
    size_t NV = 0;
    size_t NT = 0;

    for(const auto &v : mesh.tet_vertices)
        if(!v.is_removed)
            NV++;

    for(const auto &t : mesh.tets)
        if(!t.is_removed)
            NT++;

    MatrixXf TV;
    Matrix4Xi TF;

    TV.resize(3, NV);
    TF.resize(NoChange, NT);


    //Now add. Note index changes for vertices.
    std::vector<size_t> IndexTransform;
    NV = 0;

    for(const auto &v : mesh.tet_vertices)
        if(!v.is_removed){
            auto p = TV.col(NV).data();

            p[0] = (float)v[0];
            p[1] = (float)v[1];
            p[2] = (float)v[2];

            IndexTransform.push_back(NV++);
        }else
            IndexTransform.push_back(SIZE_MAX);

    NT = 0;

    for(const auto &t : mesh.tets)
        if(!t.is_removed){
            auto p = TF.col(NT++).data();

            p[0] = IndexTransform[t[0]];
            p[1] = IndexTransform[t[1]];
            p[2] = IndexTransform[t[3]];
            p[3] = IndexTransform[t[2]];
        }


    return std::make_pair(TV, TF);
}

