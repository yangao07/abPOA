cmake_minimum_required(VERSION 3.14)
project(test_parent LANGUAGES C)

add_subdirectory(${ABPOA_SOURCE_DIR} abpoa_build)

add_executable(test_parent test_parent.c)
target_link_libraries(test_parent abpoa)

# Check that abpoa's PUBLIC include directories expose include/ not src/
get_target_property(ABPOA_IFACE_DIRS abpoa INTERFACE_INCLUDE_DIRECTORIES)
string(FIND "${ABPOA_IFACE_DIRS}" "/src>" _src_pos)
if(NOT _src_pos EQUAL -1)
    message(FATAL_ERROR
        "abpoa leaks internal src/ headers to consumers via PUBLIC include: ${ABPOA_IFACE_DIRS}")
endif()

# Check that CMAKE_C_FLAGS was NOT polluted with -march=native
if(CMAKE_C_FLAGS MATCHES "-march=native")
    message(FATAL_ERROR
        "CMAKE_C_FLAGS polluted with -march=native: '${CMAKE_C_FLAGS}'")
endif()
