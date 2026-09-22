#ifndef ENERGY_H
#define ENERGY_H

#include <stddef.h>

// Use this Atom definition throughout the complete project.
typedef struct {
    char element[8];
    double x;
    double y;
    double z;
} Atom;

// Stillinger-Weber potential parameters.
typedef struct {
    double epsilon;
    double sigma;
    double cutoff_a;
    double lambda;
    double gamma;
    double cos_theta_0;
    double A;
    double B;
    double p;
    double q;
} SWParameters;

// Standard Stillinger-Weber parameters for silicon.
extern const SWParameters SILICON_SW_PARAMETERS;

// Return the number of unique atoms in a periodic cubic supercell.
//
// structure_choice:
// 1 = SC
// 2 = BCC
// 3 = FCC
// 4 = diamond cubic
size_t periodic_atom_count(
    int structure_choice,
    int nx,
    int ny,
    int nz
);

// Generate unique periodic atoms for an energy calculation.
//
// The function does not include duplicate atoms on the positive
// supercell boundaries.
size_t generate_periodic_cubic_supercell(
    Atom atoms[],
    size_t capacity,
    int structure_choice,
    const char *element,
    int nx,
    int ny,
    int nz,
    double lattice_constant
);

// Calculate the minimum-image vector from atom i to atom j.
//
// lx, ly, and lz are the complete supercell lengths.
void energy_pbc_displacement(
    const Atom *atom_i,
    const Atom *atom_j,
    double lx,
    double ly,
    double lz,
    double displacement[3]
);

// Check the minimum-image cutoff condition:
//
// 2 * cutoff < minimum box length
int minimum_image_cutoff_is_valid(
    double lx,
    double ly,
    double lz,
    double cutoff
);

// Calculate total unshifted Lennard-Jones energy.
double calculate_lj_total_energy(
    const Atom atoms[],
    size_t number_of_atoms,
    double lx,
    double ly,
    double lz,
    double epsilon,
    double sigma,
    double cutoff
);

// Calculate total Stillinger-Weber energy.
//
// The function also returns separate two-body and three-body
// contributions through the output pointers.
double calculate_sw_total_energy(
    const Atom atoms[],
    size_t number_of_atoms,
    double lx,
    double ly,
    double lz,
    const SWParameters *parameters,
    double *two_body_energy,
    double *three_body_energy
);

#endif