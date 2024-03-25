// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/FloatTetDelaunay.h>

#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

#include <iterator>
#include <algorithm>
#include <bitset>

#include <floattetwild/Predicates.hpp>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/MeshIO.hpp>

namespace floatTetWild {



uint32_t HashProbe(Mesh &mesh, Vector3 min, Vector3 max, const std::vector<Vector3> &boxpoints, const std::vector<Vector3> &voxel_points, const std::vector<double> &V_d) {
    //from http://isthe.com/chongo/src/fnv/hash_32.c

    /*
     * fnv_32_buf - perform a 32 bit Fowler/Noll/Vo hash on a buffer
     *
     * input:
     *	buf	- start of buffer to hash
     *	len	- length of buffer in octets
     *	hval	- previous hash value or 0 if first call
     *
     * returns:
     *	32 bit hash as a static hash type
     *
     * NOTE: To use the 32 bit FNV-0 historic hash, use FNV0_32_INIT as the hval
     *	 argument on the first call to either fnv_32_buf() or fnv_32_str().
     *
     * NOTE: To use the recommended 32 bit FNV-1 hash, use FNV1_32_INIT as the hval
     *	 argument on the first call to either fnv_32_buf() or fnv_32_str().
     */
    auto fnv_32_buf = [](const void *buf, size_t len, uint32_t hval = 0){
        const unsigned char *bp = (const unsigned char *)buf;	/* start of buffer */
        const unsigned char *be = bp + len;		/* beyond end of buffer */

        /*
         * FNV-1 hash each octet in the buffer
         */
        while (bp < be) {
	        /* multiply by the 32 bit FNV magic prime mod 2^32 */
	        hval *= 0x01000193;

	        /* xor the bottom with the current octet */
	        hval ^= 1 + (uint32_t)*bp++;
        }

        /* return our new hash value */
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


    hash = fnv_32_buf(min.data(), 3*sizeof(Scalar), hash);
    hash = fnv_32_buf(max.data(), 3*sizeof(Scalar), hash);
    for(const auto &e : boxpoints) hash = fnv_32_buf(e.data(), 3*sizeof(Scalar), hash);
    for(const auto &e : voxel_points) hash = fnv_32_buf(e.data(), 3*sizeof(Scalar), hash);

    hash = fnv_32_buf(V_d.data(), V_d.size()*sizeof(double), hash);

    return hash;
}








	namespace {
        void
        get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
            min = vertices.front();
            max = vertices.front();

            for (size_t j = 0; j < vertices.size(); j++) {
                for (int i = 0; i < 3; i++) {
                    min(i) = std::min(min(i), vertices[j](i));
                    max(i) = std::max(max(i), vertices[j](i));
                }
            }

//            const Scalar dis = std::max((max - min).minCoeff() * params.box_scale, params.eps_input * 2);
            const Scalar dis = std::max(params.ideal_edge_length, params.eps_input * 2);
            for (int j = 0; j < 3; j++) {
                min[j] -= dis;
                max[j] += dis;
            }

            logger().debug("min = {} {} {}", min[0], min[1], min[2]);
            logger().debug("max = {} {} {}", max[0], max[1], max[2]);
        }

        bool comp(const std::array<int, 4> &a, const std::array<int, 4> &b) {
            return std::tuple<int, int, int>(a[0], a[1], a[2]) < std::tuple<int, int, int>(b[0], b[1], b[2]);
        }

        void match_surface_fs(Mesh &mesh, const std::vector<Vector3> &input_vertices,
                              const std::vector<Vector3i> &input_faces, std::vector<bool> &is_face_inserted) {
            std::vector<std::array<int, 4>> input_fs(input_faces.size());
            for (int i = 0; i < input_faces.size(); i++) {
                input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
                std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
            }
            std::sort(input_fs.begin(), input_fs.end(), comp);

//            for(auto& f: input_fs){
//                cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<endl;
//            }
//            cout<<"/////"<<endl;

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    std::array<int, 3> f = {{t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4]}};
                    std::sort(f.begin(), f.end());
//                    cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
                    auto bounds = std::equal_range(input_fs.begin(), input_fs.end(),
                                                   std::array<int, 4>({{f[0], f[1], f[2], -1}}),
                                                   comp);
//                    bool is_matched = false;
//                    int total_ori = 0;
                    for (auto it = bounds.first; it != bounds.second; ++it) {
//                        is_matched = true;
                        int f_id = (*it)[3];
                        is_face_inserted[f_id] = true;
//                        int ori = Predicates::orient_3d(mesh.tet_vertices[t[j]].pos,
//                                                        input_vertices[input_faces[f_id][0]],
//                                                        input_vertices[input_faces[f_id][1]],
//                                                        input_vertices[input_faces[f_id][2]]);
//                        if (ori == Predicates::ORI_POSITIVE)
//                            total_ori++;
//                        else if (ori == Predicates::ORI_NEGATIVE)
//                            total_ori--;
                    }
//                    if (is_matched)
//                        t.is_surface_fs[j] = total_ori;
//                    else
//                        t.is_surface_fs[j] = NOT_SURFACE;

//                    if(is_matched)
//                        cout<<"matched: "<<total_ori<<endl;
                }
            }
        }

        void match_bbox_fs(Mesh &mesh, const Vector3 &min, const Vector3 &max) {
            auto get_bbox_fs = [&](const MeshTet &t, int j) {
                std::array<int, 6> cnts = {{0, 0, 0, 0, 0, 0}};
                for (int k = 0; k < 3; k++) {
                    Vector3 &pos = mesh.tet_vertices[t[(j + k + 1) % 4]].pos;
                    for (int n = 0; n < 3; n++) {
                        if (pos[n] == min[n])
                            cnts[n * 2]++;
                        else if (pos[n] == max[n])
                            cnts[n * 2 + 1]++;
                    }
                }
                for (int i = 0; i < cnts.size(); i++) {
                    if (cnts[i] == 3)
                        return i;
                }
                return NOT_BBOX;
            };

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    t.is_bbox_fs[j] = get_bbox_fs(t, j);
                }
            }
        }

        void
        compute_voxel_points(const Vector3 &min, const Vector3 &max, const Parameters &params, const AABBWrapper &tree,
                             std::vector<Vector3> &voxels) {
            const Vector3 diag = max - min;
            Vector3i n_voxels = (diag / (params.bbox_diag_length * params.box_scale)).cast<int>();

            for (int d = 0; d < 3; ++d)
                n_voxels(d) = std::max(n_voxels(d), 1);

            const Vector3 delta = diag.array() / n_voxels.array().cast<Scalar>();

            voxels.clear();
            voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

//            const double sq_distg = std::min(params.ideal_edge_length / 2, 10 * params.eps);
            const double sq_distg = 100 * params.eps_2;
            GEO::vec3 nearest_point;

            for (int i = 0; i <= n_voxels(0); ++i) {
                const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
                for (int j = 0; j <= n_voxels(1); ++j) {
                    const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
                    for (int k = 0; k <= n_voxels(2); ++k) {
                        const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);
                        if (tree.get_sq_dist_to_sf(Vector3(px, py, pz)) > sq_distg)
                            voxels.emplace_back(px, py, pz);
                    }
                }
            }
        }
    }

//#include <igl/unique_rows.h>
//#include <floattetwild/Predicates.hpp>
//    extern "C" floatTetWild::Scalar orient3d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc, const floatTetWild::Scalar *pd);

	void FloatTetDelaunay::tetrahedralize(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper &tree,
	        Mesh &mesh, std::vector<bool> &is_face_inserted) {
        const Parameters &params = mesh.params;
        auto &tet_vertices = mesh.tet_vertices;
        auto &tets = mesh.tets;

        is_face_inserted.resize(input_faces.size(), false);

        Vector3 min, max;
        get_bb_corners(params, input_vertices, min, max);
        mesh.params.bbox_min = min;
        mesh.params.bbox_max = max;

        std::vector<Vector3> boxpoints; //(8);


        std::vector<Vector3> voxel_points;
        compute_voxel_points(min, max, params, tree, voxel_points);

        const int n_pts = input_vertices.size() + boxpoints.size() + voxel_points.size();
        tet_vertices.resize(n_pts);
//        std::vector<double> V_d;
//        V_d.resize(n_pts * 3);

        size_t index = 0;
        int offset = 0;
        for (int i = 0; i < input_vertices.size(); i++) {
            tet_vertices[offset + i].pos = input_vertices[i];
        }
        offset += input_vertices.size();
        for (int i = 0; i < boxpoints.size(); i++) {
            tet_vertices[i + offset].pos = boxpoints[i];
        }
        offset += boxpoints.size();
        for (int i = 0; i < voxel_points.size(); i++) {
            tet_vertices[i + offset].pos = voxel_points[i];
        }

        std::vector<double> V_d;
        V_d.resize(n_pts * 3);
        for (int i = 0; i < tet_vertices.size(); i++) {
            for (int j = 0; j < 3; j++)
                V_d[i * 3 + j] = tet_vertices[i].pos[j];
        }
//cout << HashProbe(mesh, min, max, boxpoints, voxel_points, V_d) << "\n";
        GEO::Delaunay::initialize();
        GEO::Delaunay_var T = GEO::Delaunay::create(3, "BDEL");
        T->set_vertices(n_pts, V_d.data());
        //
        tets.resize(T->nb_cells());
        const auto &tet2v = T->cell_to_v();
        for (int i = 0; i < T->nb_cells(); i++) {
            for (int j = 0; j < 4; ++j) {
                const int v_id = tet2v[i * 4 + j];

                tets[i][j] = v_id;
                tet_vertices[v_id].conn_tets.push_back(i);
            }
            std::swap(tets[i][1], tets[i][3]);
        }
//cout << HashProbe(mesh, min, max, boxpoints, voxel_points, V_d) << "\n"; exit(0);
        for (int i = 0; i < mesh.tets.size(); i++) {
            auto &t = mesh.tets[i];
            if (is_inverted(mesh.tet_vertices[t[0]].pos, mesh.tet_vertices[t[1]].pos,
                            mesh.tet_vertices[t[2]].pos, mesh.tet_vertices[t[3]].pos)) {
                cout << "EXIT_INV" << endl;
                exit(0);
            }
        }



        //match faces: should be integer with sign
        //match bbox 8 facets: should be -1 and 0~5
//        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted);
        match_bbox_fs(mesh, min, max);

//        MeshIO::write_mesh("delaunay.msh", mesh);
    }
}
