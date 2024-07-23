add_library(theuerkauf_fitter TheuerkaufPeak.hpp Util.hpp)
set_target_properties(theuerkauf_fitter PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(theuerkauf_fitter ${ROOT_LIBRARIES} )
