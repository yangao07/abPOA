#!/bin/bash
# Test that libabpoa.a does not export bare klib symbols.
# All klib symbols must be prefixed with abpoa_ to avoid collisions
# with minimap2, htslib, and other klib-bundling libraries.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

LIB="lib/libabpoa.a"
if [ ! -f "$LIB" ]; then
    echo "SKIP: $LIB not found (run make first)"
    exit 1
fi

SYMBOLS=$(nm -g "$LIB" 2>/dev/null | grep " T " | awk '{print $3}')
ERRORS=0

# These bare symbols must NOT appear
BARE_SYMBOLS="kmalloc krealloc kcalloc kfree km_init km_init2 km_destroy km_stat
              kvsprintf ksprintf ksplit_core kstrtok kgetline kmemmem kstrstr kstrnstr"

for sym in $BARE_SYMBOLS; do
    if echo "$SYMBOLS" | grep -qx "$sym"; then
        echo "FAIL: bare symbol '$sym' found in $LIB"
        ERRORS=$((ERRORS + 1))
    fi
done

# These prefixed symbols MUST appear
PREFIXED_SYMBOLS="abpoa_kmalloc abpoa_krealloc abpoa_kcalloc abpoa_kfree
                   abpoa_km_init abpoa_km_init2 abpoa_km_destroy abpoa_km_stat
                   abpoa_kvsprintf abpoa_ksprintf abpoa_ksplit_core abpoa_kstrtok
                   abpoa_kgetline abpoa_kmemmem abpoa_kstrstr abpoa_kstrnstr"

for sym in $PREFIXED_SYMBOLS; do
    if ! echo "$SYMBOLS" | grep -qx "$sym"; then
        echo "FAIL: prefixed symbol '$sym' NOT found in $LIB"
        ERRORS=$((ERRORS + 1))
    fi
done

if [ $ERRORS -gt 0 ]; then
    echo "$ERRORS symbol namespacing errors"
    exit 1
fi

echo "All klib symbols properly namespaced"
exit 0
