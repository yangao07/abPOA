#!/bin/bash
# Test that SIMDMalloc/SIMDFree are portable across platforms.
# The implementation must:
# 1. Have an MSVC path with _aligned_malloc/_aligned_free
# 2. Not use a bare free() macro for SIMDFree
# 3. Compile and produce correct alignment on this platform
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

SIMD_H="src/simd_instruction.h"
ERRORS=0

if ! grep -q '_MSC_VER' "$SIMD_H"; then
    echo "FAIL: $SIMD_H has no MSVC (_MSC_VER) code path"
    ERRORS=$((ERRORS + 1))
fi

if ! grep -q '_aligned_malloc' "$SIMD_H"; then
    echo "FAIL: $SIMD_H does not use _aligned_malloc for MSVC"
    ERRORS=$((ERRORS + 1))
fi

if ! grep -q '_aligned_free' "$SIMD_H"; then
    echo "FAIL: $SIMD_H does not use _aligned_free for MSVC"
    ERRORS=$((ERRORS + 1))
fi

if grep -q '#define SIMDFree.*free' "$SIMD_H"; then
    echo "FAIL: SIMDFree is still a bare free() macro"
    ERRORS=$((ERRORS + 1))
fi

# Compile and run the alignment verification test
if [ -f tests/test_aligned_alloc.c ]; then
    if gcc -march=native tests/test_aligned_alloc.c -I src -I include \
           -o /tmp/test_aligned_alloc 2>/tmp/test_aligned_alloc.err; then
        if ! /tmp/test_aligned_alloc; then
            echo "FAIL: aligned allocation test failed at runtime"
            ERRORS=$((ERRORS + 1))
        fi
    else
        echo "FAIL: test_aligned_alloc.c does not compile"
        cat /tmp/test_aligned_alloc.err
        ERRORS=$((ERRORS + 1))
    fi
    rm -f /tmp/test_aligned_alloc /tmp/test_aligned_alloc.err
fi

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS portable allocation errors"
    exit 1
fi

echo "Portable allocation checks passed"
exit 0
