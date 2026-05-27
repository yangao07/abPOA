#!/bin/bash
# Test that CMake can build only the library (no CLI executable)
# when ABPOA_BUILD_EXE=OFF.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

BUILD_DIR="/tmp/abpoa_build_libonly_$$"
trap "rm -rf $BUILD_DIR" EXIT

ERRORS=0

# Test 1: Build with ABPOA_BUILD_EXE=OFF
cmake -B "$BUILD_DIR" -S . -DABPOA_BUILD_EXE=OFF \
      -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 || true
cmake --build "$BUILD_DIR" -j$(nproc) >/dev/null 2>&1 || true

# Library must exist
if find "$BUILD_DIR" -name 'libabpoa.a' | grep -q .; then
    : # ok
else
    echo "FAIL: libabpoa.a not found in lib-only build"
    ERRORS=$((ERRORS + 1))
fi

# Binary must NOT exist
if find "$BUILD_DIR" -name 'abpoa' -type f -executable | grep -q .; then
    echo "FAIL: abpoa binary found in lib-only build (should be disabled)"
    ERRORS=$((ERRORS + 1))
fi

# Test 2: Default top-level build should still produce the binary
BUILD_DIR2="/tmp/abpoa_build_default_$$"
trap "rm -rf $BUILD_DIR $BUILD_DIR2" EXIT

cmake -B "$BUILD_DIR2" -S . -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 || true
cmake --build "$BUILD_DIR2" -j$(nproc) >/dev/null 2>&1 || true

if ! find "$BUILD_DIR2" -name 'abpoa' -type f -executable | grep -q .; then
    echo "FAIL: abpoa binary not found in default build (should be built)"
    ERRORS=$((ERRORS + 1))
fi

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS cmake lib-only errors"
    exit 1
fi

echo "CMake library-only build works correctly"
exit 0
