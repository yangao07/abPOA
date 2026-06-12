/*
 * Test that SIMDMalloc returns properly aligned memory and SIMDFree
 * releases it without error.
 *
 * Compile:
 *   gcc -march=native tests/test_aligned_alloc.c -I src -I include -o test_aligned_alloc
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "simd_instruction.h"

int main(void) {
    int alignments[] = {16, 32, 64};
    int sizes[] = {64, 1024, 65536};

    for (int a = 0; a < 3; a++) {
        for (int s = 0; s < 3; s++) {
            void *ptr = SIMDMalloc(sizes[s], alignments[a]);
            assert(ptr != NULL);
            assert(((uintptr_t)ptr % alignments[a]) == 0);

            /* Write to verify the memory is usable */
            ((char *)ptr)[0] = 'A';
            ((char *)ptr)[sizes[s] - 1] = 'Z';

            SIMDFree(ptr);
        }
    }

    printf("PASS: SIMDMalloc/SIMDFree alignment verified\n");
    return 0;
}
