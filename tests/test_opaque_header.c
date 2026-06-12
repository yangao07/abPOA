/*
 * Test that abpoa.h can be included without SIMD flags.
 *
 * This file is compiled with only -I include and NO -march=native
 * or -msse2 flags. If abpoa.h still pulls in simd_instruction.h,
 * this will fail to compile.
 *
 * Compile (test only, no link):
 *   gcc -c -Wall -Werror -I include tests/test_opaque_header.c -o /dev/null
 */
#include <stdio.h>
#include <stdlib.h>
#include "abpoa.h"

int check_header_is_clean(void) {
    /* Verify we can declare abpoa_t pointers without SIMD knowledge */
    abpoa_t *ab = NULL;
    abpoa_para_t *abpt = NULL;
    (void)ab;
    (void)abpt;
    return 0;
}
