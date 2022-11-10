#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <tbb/parallel_for.h>
#include "rkcommon/math/vec.h"
#include <CGAL/draw_triangulation_3.h>
using namespace rkcommon::math;

// #define CGAL_LINKED_WITH_TBB
#define CGAL_USE_BASIC_VIEWER
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
    // #ifdef CGAL_LINKED_WITH_TBB


    const int NUM_INSERTED_POINTS = 50;
    
    CGAL::Random_points_in_cube_3<Point> rnd(1.);
    // Construction from a vector of 1,000,000 points
    // std::vector<Point> V;
    std::vector<std::pair<Point,int>> V;
    V.reserve(NUM_INSERTED_POINTS);
    for (int i = 0; i != NUM_INSERTED_POINTS; ++i){
        V.push_back(std::make_pair(*rnd++, i));
    }
    // for(int i = 0; i < NUM_INSERTED_POINTS; i++){
    //     Point p0 = V[i];
    //     std::cout << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
    // }
        
    // Construct the locking data-structure, using the bounding-box of the points
    Triangulation::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
    // Construct the triangulation in parallel
    Triangulation T(V.begin(), V.end(), &locking_ds);
    assert(T.is_valid());
    // CGAL::draw(T);

    Vertex_handle v;
    Cell_handle cell;
    cell = T.locate(Point(0.210878, 0.293823, 0.529068), cell);

    if(T.is_infinite(cell))
    {   
        std::cout << "Infinite cell " << std::endl;
        // seed_validity.Set(i, 0);
        // continue;
    } else{
        auto index1 = cell->vertex(0) -> info();
        // // std::cout << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
        int index2 = cell->vertex(1) ->info();
        int index3 = cell->vertex(2) ->info();
        int index4 = cell->vertex(3) ->info();
        std::cout << "debug: " << "\n";
        std::cout << index1 << " " << index2 << " " << index3 << " " << index4 << std::endl;
    }
 
   



// #endif //CGAL_LINKED_WITH_TBB
    return 0;
}