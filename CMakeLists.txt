cmake_minimum_required(VERSION 3.2)
project(abpoa LANGUAGES C VERSION 3.0.0)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.21")
    # PROJECT_IS_TOP_LEVEL is set automatically
else()
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
        set(PROJECT_IS_TOP_LEVEL TRUE)
    else()
        set(PROJECT_IS_TOP_LEVEL FALSE)
    endif()
endif()

option(ABPOA_BUILD_EXE "Build the abpoa command-line executable" ${PROJECT_IS_TOP_LEVEL})
option(ABPOA_INSTALL "Generate install targets" ${PROJECT_IS_TOP_LEVEL})

include(GNUInstallDirs)
find_package(ZLIB REQUIRED)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

# build abPOA as a static library by default
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")
if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

add_library(abpoa
    src/abpoa_align.c
    src/abpoa_graph.c
    src/abpoa_output.c
    src/abpoa_plot.c
    src/abpoa_seed.c
    src/abpoa_seq.c
    src/kalloc.c
    src/kstring.c
    src/abpoa_simd.c
    src/abpoa_align_simd.c
    src/utils.c)

target_compile_options(abpoa PRIVATE -Wall -march=native)
target_link_libraries(abpoa z pthread m)

target_include_directories(abpoa
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(ABPOA_BUILD_EXE)
    add_executable(abpoa_bin src/abpoa.c)
    target_link_libraries(abpoa_bin abpoa)
    set_target_properties(abpoa_bin PROPERTIES OUTPUT_NAME abpoa)
endif()

if(ABPOA_INSTALL)
    install(TARGETS abpoa DESTINATION ${CMAKE_INSTALL_LIBDIR})
    if(ABPOA_BUILD_EXE)
        install(TARGETS abpoa_bin DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h")

    # configure and install pkg-config file
    if(NOT BUILD_SHARED_LIBS)
        set(abPOA_ADDITIONAL_LINK_LIBRARIES "-lz -lpthread -lm")
    endif()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/abpoa.pc.in
                   ${CMAKE_CURRENT_BINARY_DIR}/abpoa-1.pc @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/abpoa-1.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif()
