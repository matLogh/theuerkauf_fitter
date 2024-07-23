
# cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
#    project(theuerkauf_fitter LANGUAGES CXX)


# add_compile_options(-std=c++17 -g -fmax-errors=5 -O3)

# # ROOT stuff, should be in the top-level CMakeLists.txt
# find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf)
# include(${ROOT_USE_FILE})
# include_directories(${ROOT_INCLUDE_DIRS})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# compile fitter lib
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/theuerkauf/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/theuerkauf/)


# add_executable(fitter_test example/test.cpp)
# # target_include_directories(fitter_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# target_link_libraries(fitter_test theuerkauf_fitter ${ROOT_LIBRARIES})
# set_target_properties(fitter_test PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
# )