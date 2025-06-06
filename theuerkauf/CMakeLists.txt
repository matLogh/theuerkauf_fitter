# add_library(theuerkauf_fitter TheuerkaufPeak.hpp Util.hpp)
# set_target_properties(theuerkauf_fitter PROPERTIES LINKER_LANGUAGE CXX)
# target_link_libraries(theuerkauf_fitter PUBLIC ${ROOT_LIBRARIES} )
# target_include_directories(theuerkauf_fitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(theuerkauf_fitter TheuerkaufPeak.cpp Util.hpp)
set_target_properties(theuerkauf_fitter PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(theuerkauf_fitter PUBLIC ${ROOT_LIBRARIES} )
target_include_directories(theuerkauf_fitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

