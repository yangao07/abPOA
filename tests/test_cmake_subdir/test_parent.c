/*
 * Minimal consumer that includes abpoa.h and uses the library API.
 * Used to verify add_subdirectory() integration works correctly.
 */
#include <stdio.h>
#include <stdlib.h>
#include "abpoa.h"

int main(void) {
    abpoa_t *ab = abpoa_init();
    if (ab == NULL) return 1;
    abpoa_free(ab);
    printf("PASS: add_subdirectory consumer works\n");
    return 0;
}
