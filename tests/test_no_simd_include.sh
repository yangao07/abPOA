#!/bin/bash
# Test that include/abpoa.h does not leak SIMD types into the public API.
# After making abpoa_simd_matrix_t opaque, the public header should:
# 1. NOT include simd_instruction.h
# 2. NOT reference SIMDi type
# 3. Compile without any -m SIMD flags
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

ERRORS=0

if grep -q '#include.*simd_instruction' include/abpoa.h; then
    echo "FAIL: include/abpoa.h still includes simd_instruction.h"
    ERRORS=$((ERRORS + 1))
fi

if grep -q 'SIMDi' include/abpoa.h; then
    echo "FAIL: include/abpoa.h still references SIMDi type"
    ERRORS=$((ERRORS + 1))
fi

# Compile a consumer file with only -I include and NO SIMD flags.
# This must succeed if the public header is clean.
if [ -f tests/test_opaque_header.c ]; then
    if ! gcc -c -Wall -Werror -I include tests/test_opaque_header.c \
         -o /tmp/test_opaque_header.o 2>/tmp/test_opaque_header.err; then
        echo "FAIL: test_opaque_header.c does not compile without SIMD flags"
        cat /tmp/test_opaque_header.err
        ERRORS=$((ERRORS + 1))
    fi
    rm -f /tmp/test_opaque_header.o /tmp/test_opaque_header.err
fi

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS opaque header errors"
    exit 1
fi

echo "Public header is SIMD-clean"
exit 0
