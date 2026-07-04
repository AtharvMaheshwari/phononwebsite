#include <emscripten.h>
#include <stdlib.h>
#include <stdio.h>
#include "spglib/spglib.h"

// We want to return rotations and translations.
// Since WASM only supports returning numbers, we will write to a buffer.
// spg_get_symmetry returns the number of operations.

EMSCRIPTEN_KEEPALIVE
int get_symmetry(
    double* lattice, // 9 doubles (3x3)
    double* position, // 3 * num_atom doubles
    int* types, // num_atom ints
    int num_atom,
    double symprec,
    int* rotations_out, // buffer for 192 * 9 ints
    double* translations_out // buffer for 192 * 3 doubles
) {
    // Spglib expects max_size = 192 for crystals with centering translations (e.g. FCC)
    int max_size = 192;

    int size_sym = spg_get_symmetry(
        (int (*)[3][3])rotations_out,
        (double (*)[3])translations_out,
        max_size,
        (const double (*)[3])lattice,
        (const double (*)[3])position,
        types,
        num_atom,
        symprec
    );
    return size_sym;
}
