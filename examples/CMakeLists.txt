add_executable(chapman chapman/chapman.cpp)

target_compile_features(chapman PUBLIC cxx_std_23)

target_link_libraries(chapman PRIVATE KinetiCpp)
target_link_libraries(chapman PRIVATE Eigen)

if (KINETICPP_SANITIZE_EXAMPLES)
    include(../cmake/Sanitizers.cmake)
    target_sanitizer_options(chapman)
else()
    
endif()