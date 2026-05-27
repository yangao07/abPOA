cmake_minimum_required(VERSION 3.10)
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

# Auto-detect architecture
if(EMSCRIPTEN OR CMAKE_SYSTEM_PROCESSOR MATCHES "^wasm")
    set(_ABPOA_WASM TRUE)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
    set(_ABPOA_X86_64 TRUE)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    set(_ABPOA_AARCH64 TRUE)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^armv7")
    set(_ABPOA_ARM32 TRUE)
endif()

if(_ABPOA_X86_64)
    option(ABPOA_SIMD_DISPATCH "Build with runtime SIMD dispatch (x86_64)" ON)
else()
    set(ABPOA_SIMD_DISPATCH OFF)
endif()

include(GNUInstallDirs)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

# build abPOA as a static library by default
if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS OFF)
endif()
if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

# Base sources (non-SIMD)
set(ABPOA_BASE_SOURCES
    src/abpoa_align.c
    src/abpoa_graph.c
    src/abpoa_output.c
    src/abpoa_plot.c
    src/abpoa_seed.c
    src/abpoa_seq.c
    src/kalloc.c
    src/kstring.c
    src/abpoa_simd.c
    src/utils.c)

add_library(abpoa ${ABPOA_BASE_SOURCES})
target_compile_options(abpoa PRIVATE -Wall -Wno-unused-function)
target_link_libraries(abpoa ZLIB::ZLIB Threads::Threads m)

target_include_directories(abpoa
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)

# SIMDE definitions only needed when NOT using x86_64 multi-ISA dispatch
# (dispatch mode compiles SIMD sources as separate OBJECT libraries
# that carry their own definitions; base sources use native <immintrin.h>)
if(NOT (ABPOA_SIMD_DISPATCH AND _ABPOA_X86_64))
    target_compile_definitions(abpoa PRIVATE USE_SIMDE SIMDE_ENABLE_NATIVE_ALIASES)
endif()

# --- SIMD compilation strategy ---
if(ABPOA_SIMD_DISPATCH AND _ABPOA_X86_64)
    # Propagate PIC to OBJECT libraries so shared library builds work
    get_target_property(_abpoa_pic abpoa POSITION_INDEPENDENT_CODE)

    # Multi-ISA dispatch: compile align source 4x (once per ISA)
    set(_ABPOA_ISA_NAMES  sse2    sse41    avx2    avx512bw)
    set(_ABPOA_ISA_FLAGS  "-msse2" "-msse4.1" "-mavx2" "-mavx512bw")
    set(_ABPOA_ISA_EXTRAS "-U__SSE4_1__" "" "" "")

    list(LENGTH _ABPOA_ISA_NAMES  _isa_count)
    list(LENGTH _ABPOA_ISA_FLAGS  _n_flags)
    list(LENGTH _ABPOA_ISA_EXTRAS _n_extras)
    if(NOT _isa_count EQUAL _n_flags OR NOT _isa_count EQUAL _n_extras)
        message(FATAL_ERROR "ISA name/flag/extra lists have different lengths")
    endif()
    math(EXPR _isa_last "${_isa_count} - 1")

    foreach(_idx RANGE 0 ${_isa_last})
        list(GET _ABPOA_ISA_NAMES  ${_idx} _isa_name)
        list(GET _ABPOA_ISA_FLAGS  ${_idx} _isa_flag)
        list(GET _ABPOA_ISA_EXTRAS ${_idx} _isa_extra)

        add_library(abpoa_align_${_isa_name} OBJECT src/abpoa_align_simd.c)
        set_target_properties(abpoa_align_${_isa_name} PROPERTIES
            POSITION_INDEPENDENT_CODE "${_abpoa_pic}")
        target_compile_options(abpoa_align_${_isa_name} PRIVATE
            -Wall -Wno-unused-function ${_isa_flag})
        if(_isa_extra)
            target_compile_options(abpoa_align_${_isa_name} PRIVATE ${_isa_extra})
        endif()
        target_compile_definitions(abpoa_align_${_isa_name} PRIVATE
            ABPOA_SIMD_DISPATCH USE_SIMDE SIMDE_ENABLE_NATIVE_ALIASES)
        target_include_directories(abpoa_align_${_isa_name} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/src)

        target_sources(abpoa PRIVATE
            $<TARGET_OBJECTS:abpoa_align_${_isa_name}>)
    endforeach()

    # Dispatch source compiled once — contains only the runtime
    # CPUID dispatcher, no ISA-specific instructions
    add_library(abpoa_dispatch OBJECT src/abpoa_dispatch_simd.c)
    set_target_properties(abpoa_dispatch PROPERTIES
        POSITION_INDEPENDENT_CODE "${_abpoa_pic}")
    target_compile_options(abpoa_dispatch PRIVATE -Wall -Wno-unused-function)
    target_compile_definitions(abpoa_dispatch PRIVATE
        ABPOA_SIMD_DISPATCH USE_SIMDE SIMDE_ENABLE_NATIVE_ALIASES)
    target_include_directories(abpoa_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_sources(abpoa PRIVATE $<TARGET_OBJECTS:abpoa_dispatch>)

elseif(_ABPOA_AARCH64)
    # ARM64: single compile with SIMDE translating AVX2 to NEON
    target_sources(abpoa PRIVATE src/abpoa_align_simd.c)
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        target_compile_options(abpoa PRIVATE -mcpu=native)
    else()
        target_compile_options(abpoa PRIVATE -march=armv8-a+simd)
    endif()
    target_compile_definitions(abpoa PRIVATE __AVX2__)

elseif(_ABPOA_ARM32)
    # ARMv7 32-bit: single compile with NEON
    target_sources(abpoa PRIVATE src/abpoa_align_simd.c)
    target_compile_options(abpoa PRIVATE -march=armv7-a -mfpu=neon)
    target_compile_definitions(abpoa PRIVATE __AVX2__)

elseif(_ABPOA_WASM)
    # WASM SIMD128: 128-bit vectors, maps to SSE2 path through SIMDE
    # (USE_SIMDE and SIMDE_ENABLE_NATIVE_ALIASES set on base target above)
    target_sources(abpoa PRIVATE src/abpoa_align_simd.c)
    target_compile_options(abpoa PRIVATE -msimd128)

else()
    # Fallback: single compile with -march=native
    target_sources(abpoa PRIVATE src/abpoa_align_simd.c)
    target_compile_options(abpoa PRIVATE -march=native)
endif()

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
