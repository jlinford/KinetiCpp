include(../cmake/Sanitizers.cmake)

add_executable(chapman chapman/chapman.cpp)

target_compile_features(chapman PUBLIC cxx_std_23)

target_link_libraries(chapman PRIVATE KinetiCpp)
target_link_libraries(chapman PRIVATE Eigen)

target_sanitizer_options(chapman)