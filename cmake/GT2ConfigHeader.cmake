include(CheckCXXSourceCompiles)

SET(CMAKE_REQUIRED_FLAGS ${CXX_FLAGS})
check_cxx_source_compiles(
"
#include <memory>
int main(){auto a = std::make_unique<int>();}
" GENETRAIL2_HAS_MAKE_UNIQUE)

configure_file(cmake/config.h.in "${PROJECT_BINARY_DIR}/libraries/genetrail2/core/config.h")
