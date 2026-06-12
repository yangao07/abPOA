/*
 * Test that abPOA's internal allocator uses namespaced symbols.
 *
 * This program defines its own kmalloc/kfree. If abPOA still exports
 * bare kmalloc/kfree, the linker will either error (multiple definitions)
 * or silently pick one — both prove the collision problem.
 *
 * After namespacing, abPOA uses abpoa_kmalloc internally, so our
 * kmalloc/kfree are never called and the flag stays 0.
 *
 * Compile:
 *   gcc tests/test_dual_link.c -Iinclude -Llib -labpoa -lz -lm -lpthread \
 *       -o tests/test_dual_link
 * Run:
 *   ./tests/test_dual_link
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "abpoa.h"

static int foreign_kmalloc_called = 0;

void *kmalloc(void *km, size_t size) {
    foreign_kmalloc_called = 1;
    (void)km;
    return malloc(size);
}

void kfree(void *km, void *ptr) {
    foreign_kmalloc_called = 1;
    (void)km;
    free(ptr);
}

int main(void) {
    abpoa_t *ab = abpoa_init();
    assert(ab != NULL);
    abpoa_free(ab);

    assert(foreign_kmalloc_called == 0 &&
           "abPOA called bare kmalloc/kfree instead of abpoa_kmalloc/abpoa_kfree");

    printf("PASS: abPOA uses namespaced klib symbols, no collision\n");
    return 0;
}
