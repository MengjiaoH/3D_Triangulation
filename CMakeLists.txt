cmake_minimum_required(VERSION 3.16)
project(3D_Triangulation)
# set(CMAKE_C_COMPILER /usr/local/bin/gcc-11)
# set(CMAKE_CXX_COMPILER /usr/local/bin/g++-11)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# VTKm
# set(VTKm_DIR "/Users/hanmj/Documents/externals/vtk-m-v1.7.1/install/lib/cmake/vtkm-1.7")
set(VTKm_DIR "/home/mengjiao/Desktop/externals/vtk-m-v1.7.1/install/lib/cmake/vtkm-1.7")
find_package(VTKm REQUIRED)

## VTK
# set(VTK_DIR "/Users/hanmj/Documents/externals/VTK-9.1.0/install/lib/cmake/vtk-9.1")
# set(VTK_DIR "/home/mengjiao/externals/VTK-9.1.0/install/lib/cmake/vtk-9.1")
# find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})

## TBB
set(TBB_DIR "/home/mengjiao/Desktop/externals/oneTBB/install/lib/cmake/TBB")
# set(TBB_DIR "/Users/hanmj/Documents/externals/oneTBB/install/lib/cmake/TBB")

## rkcommon
# set(rkcommon_DIR "/Users/hanmj/Documents/externals/rkcommon/install/lib/cmake/rkcommon-1.10.0")
set(rkcommon_DIR "/home/mengjiao/Desktop/externals/rkcommon/install/lib/cmake/rkcommon-1.10.0")
find_package(rkcommon REQUIRED)

find_package(CGAL)
include(CGAL_TBB_support) 

# add_subdirectory(FTLE)
# add_subdirectory(AdaptiveSampling)

add_subdirectory("${CMAKE_SOURCE_DIR}/helpers")
include_directories("${CMAKE_SOURCE_DIR}/helpers")

add_executable(2D_Triangulation 2D_Triangulation.cpp)
target_include_directories(2D_Triangulation PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(2D_Triangulation PRIVATE rkcommon::rkcommon TBB::tbb helpers vtkm_cont vtkm_io ${VTKm_LIBRARIES} CGAL::CGAL CGAL::TBB_support)
target_compile_options(2D_Triangulation PRIVATE ${VTKm_COMPILE_OPTIONS})
set_property(TARGET 2D_Triangulation PROPERTY CXX_STANDARD 17)

add_executable(3D_Triangulation 3D_Triangulation.cpp)
target_include_directories(3D_Triangulation PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(3D_Triangulation PRIVATE rkcommon::rkcommon TBB::tbb helpers vtkm_cont vtkm_io ${VTKm_LIBRARIES} CGAL::CGAL CGAL::TBB_support)
target_compile_options(3D_Triangulation PRIVATE ${VTKm_COMPILE_OPTIONS})
set_property(TARGET 3D_Triangulation PROPERTY CXX_STANDARD 17)

# add_executable(double_gyre_ftle simu_dg_ftle.cpp)
# target_include_directories(double_gyre_ftle PRIVATE ${VTKm_INCLUDE_DIRS})
# target_link_libraries(double_gyre_ftle PRIVATE vtkm_cont vtkm_io ${VTKm_LIBRARIES} ${VTK_LIBRARIES} rkcommon::rkcommon TBB::tbb)
# target_compile_options(double_gyre_ftle PRIVATE ${VTKm_COMPILE_OPTIONS})
# set_property(TARGET double_gyre_ftle PROPERTY CXX_STANDARD 17)
