#pragma once
#include <floattetwild/Mesh.hpp>


//Like the exe but, you know, a library call. Takes and returns primitive data.
std::pair<Eigen::MatrixXf, Eigen::Matrix4Xi> fTetWildRun(const Eigen::MatrixXf &Vertices, const Eigen::Matrix3Xi &Faces, floatTetWild::Scalar ideal_edge_length_rel = 0.05f, std::string OutputName = "", bool Quiet = true, unsigned int max_threads = UINT_MAX);

