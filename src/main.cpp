// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <CLI/CLI.hpp>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
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

using namespace floatTetWild;
using namespace Eigen;

/*
Mesh fTetWildRun(){
    Mesh        mesh;
    Parameters& params = mesh.params;

    return mesh;
}*/

int main(int argc, char** argv){

#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    //Note: this is necessary and results in printing, at least in debug mode :(
    //It's inconsistent.
    GEO::initialize();

    bool skip_simplify = false;

    Mesh        mesh;
    Parameters& params = mesh.params;
    params.log_level = 6;
    params.is_quiet = true;


    CLI::App command_line {"float-tetwild"};
    command_line
      .add_option("-i,--input",
                  params.input_path,
                  "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
      ->check(CLI::ExistingFile);
    command_line.add_option("-o,--output",
                            params.output_path,
                            "Output tetmesh OUTPUT in .msh format. (string, optional, default: "
                            "input_file+postfix+'.msh')");

    command_line.add_option(
      "-l,--lr",
      params.ideal_edge_length_rel,
      "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");
    command_line.add_option("-e,--epsr",
                            params.eps_rel,
                            "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)");


#ifdef NEW_ENVELOPE
    std::string epsr_tags;
    command_line
      .add_option("--epsr-tags", epsr_tags, "List of envelope size for each input faces.")
      ->check(CLI::ExistingFile);
#endif

    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
#ifdef FLOAT_TETWILD_USE_TBB
    command_line.add_option("--max-threads", max_threads, "Maximum number of threads used");
#endif

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }

#ifdef FLOAT_TETWILD_USE_TBB
    const size_t MB          = 1024 * 1024;
    const size_t stack_size  = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(max_threads, num_threads);
    params.num_threads       = num_threads;

    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif


    Logger::init(!params.is_quiet, params.log_path);
    params.log_level = std::max(0, std::min(6, params.log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));


     if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    std::string output_mesh_name = params.output_path;
    if (params.output_path.size() > 3 &&
        params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) ==
          "msh")
        output_mesh_name = params.output_path;
    else if (params.output_path.size() > 4 &&
             params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) ==
               "mesh")
        output_mesh_name = params.output_path;
    else
        output_mesh_name = params.output_path + "_" + params.postfix + ".msh";


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
    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;


    /// set envelope
    GEO::Mesh                sf_mesh;
    json                     tree_with_ids;
    std::vector<std::string> meshes;
    {
#ifdef NEW_ENVELOPE
        if (!MeshIO::load_mesh(params.input_path,
                               input_vertices,
                               input_faces,
                               sf_mesh,
                               input_tags,
                               params.input_epsr_tags)) {
#else
        if (!MeshIO::load_mesh(
              params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
#endif
            logger().error("Unable to load mesh at {}", params.input_path);
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }
        else if (input_vertices.empty() || input_faces.empty()) {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }

        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
    }
    AABBWrapper tree(sf_mesh);
    if (!params.init(tree.get_sf_diag())) {
        return EXIT_FAILURE;
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty())
        tree.init_sf_tree(
          input_vertices, input_faces, params.input_epsr_tags, params.bbox_diag_length);
    else
        tree.init_sf_tree(input_vertices, input_faces, params.eps);
#endif


    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);

    if (params.log_level <= 1)
        output_component(input_vertices, input_faces, input_tags);

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




    MeshIO::write_mesh(output_mesh_name, mesh, false);


//return EXIT_SUCCESS;
//Note: this computes/writes boundary surface of tetrahedral mesh
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    if (params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    igl::write_triangle_mesh(params.output_path + "_" + params.postfix + "_sf.obj", V_sf, F_sf);


    return EXIT_SUCCESS;
}

