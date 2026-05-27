#!/bin/bash
# Test that the CMake build produces a portable binary with runtime
# SIMD dispatch, matching the Makefile build's multi-ISA compilation.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

BUILD_DIR="/tmp/abpoa_cmake_dispatch_$$"
trap "rm -rf $BUILD_DIR" EXIT

ERRORS=0

# Build via CMake
if ! cmake -B "$BUILD_DIR" -S . -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1; then
    echo "FAIL: cmake configure failed"
    ERRORS=$((ERRORS + 1))
fi

if [ $ERRORS -eq 0 ]; then
    if ! cmake --build "$BUILD_DIR" -j$(nproc) >/dev/null 2>&1; then
        echo "FAIL: cmake build failed"
        ERRORS=$((ERRORS + 1))
    fi
fi

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS cmake dispatch errors"
    exit 1
fi

# Find the library
CMAKE_LIB=$(find "$BUILD_DIR" -name 'libabpoa.a' | head -1)
if [ -z "$CMAKE_LIB" ]; then
    echo "FAIL: libabpoa.a not found"
    exit 1
fi

SYMBOLS=$(nm -g "$CMAKE_LIB" 2>/dev/null | grep " T " | awk '{print $3}')

# Check for per-ISA dispatch symbols
for sym in simd_sse2_abpoa_align_sequence_to_subgraph \
           simd_sse41_abpoa_align_sequence_to_subgraph \
           simd_avx2_abpoa_align_sequence_to_subgraph \
           simd_avx512_abpoa_align_sequence_to_subgraph \
           simd_abpoa_align_sequence_to_subgraph \
           simd_abpoa_align_sequence_to_graph \
           x86_simd; do
    if ! echo "$SYMBOLS" | grep -qx "$sym"; then
        echo "FAIL: dispatch symbol '$sym' not found in CMake-built library"
        ERRORS=$((ERRORS + 1))
    fi
done

# Compare output with Makefile reference
CMAKE_BIN=$(find "$BUILD_DIR" -name 'abpoa' -type f -executable | head -1)
if [ -z "$CMAKE_BIN" ]; then
    echo "FAIL: abpoa binary not found in CMake build"
    ERRORS=$((ERRORS + 1))
else
    if ! "$CMAKE_BIN" test_data/seq.fa 2>/dev/null | diff -q - tests/ref_consensus.txt >/dev/null 2>&1; then
        echo "FAIL: CMake binary output differs from Makefile reference (consensus)"
        ERRORS=$((ERRORS + 1))
    fi

    if ! "$CMAKE_BIN" test_data/seq.fa -a1 2>/dev/null | diff -q - tests/ref_msa.txt >/dev/null 2>&1; then
        echo "FAIL: CMake binary output differs from Makefile reference (MSA)"
        ERRORS=$((ERRORS + 1))
    fi

    if ! "$CMAKE_BIN" test_data/heter.fa -d2 2>/dev/null | diff -q - tests/ref_heter.txt >/dev/null 2>&1; then
        echo "FAIL: CMake binary output differs from Makefile reference (heter)"
        ERRORS=$((ERRORS + 1))
    fi
fi

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS cmake dispatch errors"
    exit 1
fi

echo "CMake SIMD dispatch build works correctly"
exit 0
