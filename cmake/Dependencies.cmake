include(cmake/CPM.cmake)

CPMAddPackage(NAME Eigen VERSION 3.4.0 DOWNLOAD_ONLY YES
              URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)
if(Eigen_ADDED)
  add_library(Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

