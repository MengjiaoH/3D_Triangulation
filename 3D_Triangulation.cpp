#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <tbb/parallel_for.h>
#include "rkcommon/math/vec.h"
#include "writer.h"
#include "place_seeds.h"
#include "bc_tet.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace rkcommon::math;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb; 
// Delaunay T3
// typedef CGAL::Triangulation_data_structure_3<
//         CGAL::Triangulation_vertex_base_3<K>,
//         CGAL::Delaunay_triangulation_cell_base_3<K>,
//         CGAL::Parallel_tag>                               Tds;
typedef CGAL::Triangulation_data_structure_3<Vb, CGAL::Triangulation_cell_base_3<K>, CGAL::Parallel_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>              Triangulation;
typedef Triangulation::Point                                Point;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef CGAL::Creator_uniform_3<double,K::Point_3>          Creator;




int main(int argc, char **argv)
{
    vec3i dims = vec3i(128, 128, 128);
    vec2f x_range = vec2f(0, 2 * M_PI);
    vec2f y_range = vec2f(0, 2 * M_PI);
    vec2f z_range = vec2f(0, 2 * M_PI);

    int num_new_seeds = 1000;
    int num_basis_seeds = dims.x * dims.y * dims.z;


    //** Load basic flow maps 
    auto start_loading = std::chrono::high_resolution_clock::now();
    std::string file_dir = "/home/mengjiao/Desktop/3D_Triangulation/datasets/abc/128x128x128/";
    std::vector<std::filesystem::path> filenames;
    for (const auto& entry : std::filesystem::directory_iterator{file_dir}) {
        if (entry.is_regular_file() && (entry.path().extension() == ".txt")) {
            filenames.push_back(entry.path().filename());
        }
    }
    std::sort(filenames.begin(), filenames.end(), [](const auto& lhs, const auto& rhs) { return lhs.string() < rhs.string();});
    // for(int i = 0; i < filenames.size(); ++i){
    //     std::cout << filenames[i] << std::endl;
    // }

    int num_fm = filenames.size() - 1;
    std::vector<std::vector<vec3f>> flow_maps(filenames.size(), std::vector<vec3f>(num_basis_seeds, vec3f(0.f)));

    for(int n = 0; n < filenames.size(); ++n){
        std::string fdir = file_dir + filenames[n].string();
        flow_maps[n] = read_vec3_from_txt(fdir, false);
    }
    // for(int i = 0; i < flow_maps[0].size(); ++i){
    //     std::cout << flow_maps[0][i] << std::endl;
    // }
    auto stop_loading = std::chrono::high_resolution_clock::now();
    std::cout << "Loading basis flow maps took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop_loading - start_loading).count()
              << " milliseconds\n";

    //** Construct Triangulation
    auto start_construction = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<Point,int>> current_basis;
    for(int s = 0; s < flow_maps[0].size(); ++s){
        current_basis.push_back(std::make_pair(Point(flow_maps[0][s].x, flow_maps[0][s].y, flow_maps[0][s].z), s));
    }

    // Construct the locking data-structure, using the bounding-box of the points
    Triangulation::Lock_data_structure locking_ds(CGAL::Bbox_3(x_range.x, y_range.x, z_range.x, x_range.y, y_range.y, z_range.y), 50);
    // Construct the triangulation in parallel
    Triangulation T(current_basis.begin(), current_basis.end(), &locking_ds);
    assert(T.is_valid());

    auto stop_construction = std::chrono::high_resolution_clock::now();
    std::cout << "Triangulation took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop_construction - start_construction).count()
              << " milliseconds\n";
   

    //** Generate new seeds 
    std::vector<vec3f> seeds = place_random_seeds_3d(x_range, y_range, z_range, num_new_seeds);
    // std::vector<vec3f> seeds;
    // seeds.push_back(vec3f(1.2f, 2.0f, 2.5f));
    std::vector<std::vector<vec3f>> new_trajs(num_new_seeds, std::vector<vec3f>(num_fm, vec3f(0)));
    
    auto start_interpolation = std::chrono::high_resolution_clock::now();
    //** Find interpolation neighborhood
    // for(int i = 0; i < num_new_seeds; ++i){
    //     vec3f current_seed = seeds[i];
    //     Vertex_handle v;
    //     Cell_handle cell;
    //     cell = T.locate(Point(current_seed.x, current_seed.y, current_seed.z), cell);
    //     if(T.is_infinite(cell)){   
    //         std::cout << "Infinite cell " << std::endl;
    //         continue;
    //     } else{
    //         int index1 = cell->vertex(0) -> info();
    //         int index2 = cell->vertex(1) ->info();
    //         int index3 = cell->vertex(2) ->info();
    //         int index4 = cell->vertex(3) ->info();
    //         // std::cout << "index: " << i << "\n";
    //         // std::cout << index1 << " " << index2 << " " << index3 << " " << index4 << std::endl;
    //         // std::cout << "interpolations: " << "\n";
    //         // std::cout << flow_maps[0][index1] << "\n";
    //         // std::cout << flow_maps[0][index2] << "\n";
    //         // std::cout << flow_maps[0][index3] << "\n";
    //         // std::cout << flow_maps[0][index4] << "\n";
    //         vec4f BC;
    //         vec3f p = current_seed;
    //         vec3f a = flow_maps[0][index1];
    //         vec3f b = flow_maps[0][index2];
    //         vec3f c = flow_maps[0][index3];
    //         vec3f d = flow_maps[0][index4];

    //         bary_tet(BC, p, a, b, c, d); 

    //         // start interpolate through flow maps
    //         for(int f = 1; f < flow_maps.size(); ++f){
    //             std::vector<vec3f> current_fm = flow_maps[f];
    //             vec3f interp;
    //             interp.x = BC.x * current_fm[index1].x + BC.y * current_fm[index2].x + BC.z * current_fm[index3].x + BC.w * current_fm[index4].x;
    //             interp.y = BC.x * current_fm[index1].y + BC.y * current_fm[index2].y + BC.z * current_fm[index3].y + BC.w * current_fm[index4].y;
    //             interp.z = BC.x * current_fm[index1].z + BC.y * current_fm[index2].z + BC.z * current_fm[index3].z + BC.w * current_fm[index4].z;
    //             // std::cout << "interpolation: " << interp << "\n";
    //             new_trajs[i][f-1] = interp;
    //         }// end of interpolation
    //     } // end check if cell is in T
        
    // } // end for loop of seeds

    tbb::parallel_for( tbb::blocked_range<int>(0, num_new_seeds),[&](tbb::blocked_range<int> r){
        for (int i=r.begin(); i<r.end(); ++i)
        {
            vec3f current_seed = seeds[i];
            Vertex_handle v;
            Cell_handle cell;
            cell = T.locate(Point(current_seed.x, current_seed.y, current_seed.z), cell);
            if(T.is_infinite(cell)){   
                std::cout << "Infinite cell " << std::endl;
                // continue;
            } else{
                int index1 = cell->vertex(0) -> info();
                int index2 = cell->vertex(1) ->info();
                int index3 = cell->vertex(2) ->info();
                int index4 = cell->vertex(3) ->info();
                // std::cout << "index: " << i << "\n";
                // std::cout << index1 << " " << index2 << " " << index3 << " " << index4 << std::endl;
                // std::cout << "interpolations: " << "\n";
                // std::cout << flow_maps[0][index1] << "\n";
                // std::cout << flow_maps[0][index2] << "\n";
                // std::cout << flow_maps[0][index3] << "\n";
                // std::cout << flow_maps[0][index4] << "\n";
                vec4f BC;
                vec3f p = current_seed;
                vec3f a = flow_maps[0][index1];
                vec3f b = flow_maps[0][index2];
                vec3f c = flow_maps[0][index3];
                vec3f d = flow_maps[0][index4];

                bary_tet(BC, p, a, b, c, d); 

                // 
                for(int j = 1; j < flow_maps.size(); ++j){
                    std::vector<vec3f> current_fm = flow_maps[j];
                    vec3f interp;
                    interp.x = BC.x * current_fm[index1].x + BC.y * current_fm[index2].x + BC.z * current_fm[index3].x + BC.w * current_fm[index4].x;
                    interp.y = BC.x * current_fm[index1].y + BC.y * current_fm[index2].y + BC.z * current_fm[index3].y + BC.w * current_fm[index4].y;
                    interp.z = BC.x * current_fm[index1].z + BC.y * current_fm[index2].z + BC.z * current_fm[index3].z + BC.w * current_fm[index4].z;
                    // std::cout << "interpolation: " << interp << "\n";
                    new_trajs[i][j-1] = interp;

                };

            } // end check if cell is in T
        }
    });
        
        
    
    auto stop_interpolation = std::chrono::high_resolution_clock::now();
    std::cout << "Interpolation took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop_interpolation - start_interpolation).count()
              << " milliseconds\n";

    // save the new trajectories 


    return 0;
}