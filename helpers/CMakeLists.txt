
add_library(helpers place_seeds.cpp 
                    writer.cpp
					sobol.cpp
					bc_tet.cpp)

target_link_libraries(helpers PUBLIC rkcommon::rkcommon)
set_target_properties(helpers PROPERTIES
	CXX_STANDARD 17
	CXX_STANDARD_REQUIRED ON)