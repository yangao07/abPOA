#!/bin/bash
# Test that abPOA can be used via add_subdirectory() without polluting
# the parent project's compile flags or install rules.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

BUILD_DIR="/tmp/abpoa_subdir_build_$$"
trap "rm -rf $BUILD_DIR" EXIT

ERRORS=0

# Configure the parent project that add_subdirectory()'s abPOA
if ! cmake -B "$BUILD_DIR" -S tests/test_cmake_subdir \
           -DABPOA_SOURCE_DIR="$PROJECT_DIR" \
           -DCMAKE_BUILD_TYPE=Release 2>/tmp/cmake_subdir_cfg.log; then
    echo "FAIL: cmake configure failed for add_subdirectory test"
    cat /tmp/cmake_subdir_cfg.log
    ERRORS=$((ERRORS + 1))
fi

# Build
if [ $ERRORS -eq 0 ]; then
    if ! cmake --build "$BUILD_DIR" -j$(nproc) 2>/tmp/cmake_subdir_build.log; then
        echo "FAIL: cmake build failed for add_subdirectory test"
        cat /tmp/cmake_subdir_build.log
        ERRORS=$((ERRORS + 1))
    fi
fi

# Run the consumer binary
if [ $ERRORS -eq 0 ]; then
    if ! "$BUILD_DIR/test_parent" >/dev/null 2>&1; then
        echo "FAIL: test_parent binary failed"
        ERRORS=$((ERRORS + 1))
    fi
fi

rm -f /tmp/cmake_subdir_cfg.log /tmp/cmake_subdir_build.log

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS add_subdirectory errors"
    exit 1
fi

echo "add_subdirectory integration works correctly"
exit 0
