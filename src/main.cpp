// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <CLI/CLI.hpp>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/MeshIO.hpp>

#include <geogram/mesh/mesh.h>
#include "fTetWildRun.h"

using namespace floatTetWild;
using namespace Eigen;



int main(int argc, char** argv){
    GEO::initialize();
    CLI::App command_line {"float-tetwild"};

    std::string input_path;
    floatTetWild::Scalar ideal_edge_length_rel;

    command_line.add_option("-i,--input", input_path, "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->check(CLI::ExistingFile);
    command_line.add_option("-l,--lr", ideal_edge_length_rel, "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");


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


    //Load input mesh.
    GEO::Mesh             sf_mesh;
    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;

    if (!MeshIO::load_mesh(input_path, input_vertices, input_faces, sf_mesh, input_tags))
        return EXIT_FAILURE;
    else if (input_vertices.empty() || input_faces.empty())
        return EXIT_FAILURE;


    //Convert above loaded stuff to minimalistic_mesh.
    Eigen::MatrixXf Vertices;
    Eigen::Matrix3Xi Faces;

    Vertices.resize(3, input_vertices.size());
    Faces.resize(NoChange, input_faces.size());

    for(size_t i = 0; i != input_vertices.size(); i++) Vertices.block(0, i, 3, 1) = input_vertices[i].cast<float>();
    for(size_t i = 0; i !=    input_faces.size(); i++)    Faces.block(0, i, 3, 1) = input_faces[i];


    //Tetrahedralize.
    auto mesh = fTetWildRun(Vertices, Faces, ideal_edge_length_rel, "output", true, max_threads);


    return EXIT_SUCCESS;
}




