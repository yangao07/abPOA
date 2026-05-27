#!/bin/bash
# Test that abPOA's klib symbols are properly namespaced by compiling
# a program that defines its own bare kmalloc/kfree and linking against
# libabpoa.a. If abPOA still exports bare symbols, the linker will
# either error or our test's assertion will catch the collision.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

LIB="lib/libabpoa.a"
if [ ! -f "$LIB" ]; then
    echo "SKIP: $LIB not found (run make first)"
    exit 1
fi

BUILD_DIR="/tmp/abpoa_dual_link_$$"
trap "rm -rf $BUILD_DIR" EXIT
mkdir -p "$BUILD_DIR"

if ! gcc tests/test_dual_link.c -Iinclude -Llib -labpoa -lz -lm -lpthread \
       -o "$BUILD_DIR/test_dual_link" 2>"$BUILD_DIR/compile.log"; then
    echo "FAIL: test_dual_link.c failed to compile"
    cat "$BUILD_DIR/compile.log"
    exit 1
fi

if ! "$BUILD_DIR/test_dual_link" 2>"$BUILD_DIR/run.log"; then
    echo "FAIL: test_dual_link assertion failed"
    cat "$BUILD_DIR/run.log"
    exit 1
fi

echo "Dual-link symbol collision test passed"
exit 0
