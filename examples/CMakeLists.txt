add_executable(chapman_box chapman/box_model.cpp)
target_compile_features(chapman_box PUBLIC cxx_std_23)
target_link_libraries(chapman_box PRIVATE KinetiCpp)
target_link_libraries(chapman_box PRIVATE Eigen)

add_executable(chapman_bench chapman/benchmark.cpp)
target_compile_features(chapman_bench PUBLIC cxx_std_23)
target_link_libraries(chapman_bench PRIVATE KinetiCpp)
target_link_libraries(chapman_bench PRIVATE Eigen)


if (KINETICPP_SANITIZE_EXAMPLES)
    include(../cmake/Sanitizers.cmake)
    target_sanitizer_options(chapman_box)
    target_sanitizer_options(chapman_bench)
else()
    
endif()