INCLUDE(CheckCXXSourceCompiles) 

# Check for constexpr std::abs()
check_cxx_source_compiles("
#include <cstdlib>
int main() {
    constexpr bool nonzero = std::abs(1.0) > 0;
    return 0;
}
"
HAVE_CONSTEXPR_STD_ABS)
if (NOT HAVE_CONSTEXPR_STD_ABS)
    message(FATAL_ERROR "C++ compiler '${CMAKE_CXX_COMPILER}' does not support constexpr std::abs()")
endif()
