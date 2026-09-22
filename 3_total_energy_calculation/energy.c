#include "energy.h"

#include <math.h>
#include <stdio.h>

// Return the number of basis atoms in one conventional cubic cell.
static int conventional_basis_size(int structure_choice)
{
    switch (structure_choice) {
        case 1:
            return 1;

        case 2:
            return 2;

        case 3:
            return 4;

        case 4:
            return 8;

        default:
            return 0;
    }
}

// Return the number of unique periodic atoms.
size_t periodic_atom_count(
    int structure_choice,
    int nx,
    int ny,
    int nz
)
{
    int basis_size =
        conventional_basis_size(structure_choice);

    if (basis_size == 0 ||
        nx <= 0 ||
        ny <= 0 ||
        nz <= 0) {

        return 0;
    }

    return (size_t) basis_size
         * (size_t) nx
         * (size_t) ny
         * (size_t) nz;
}

// Generate a unique periodic cubic supercell.
size_t generate_periodic_cubic_supercell(
    Atom atoms[],
    size_t capacity,
    int structure_choice,
    const char *element,
    int nx,
    int ny,
    int nz,
    double lattice_constant
)
{
    // Simple-cubic basis.
    static const double sc_basis[1][3] = {
        {0.0, 0.0, 0.0}
    };

    // Body-centred-cubic basis.
    static const double bcc_basis[2][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    // Face-centred-cubic basis.
    static const double fcc_basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0}
    };

    // Diamond-cubic conventional-cell basis.
    static const double diamond_basis[8][3] = {
        {0.0,  0.0,  0.0},
        {0.0,  0.5,  0.5},
        {0.5,  0.0,  0.5},
        {0.5,  0.5,  0.0},
        {0.25, 0.25, 0.25},
        {0.25, 0.75, 0.75},
        {0.75, 0.25, 0.75},
        {0.75, 0.75, 0.25}
    };

    const double (*basis)[3] = NULL;
    int basis_size = 0;

    if (atoms == NULL ||
        element == NULL ||
        nx <= 0 ||
        ny <= 0 ||
        nz <= 0 ||
        lattice_constant <= 0.0) {

        return 0;
    }

    switch (structure_choice) {
        case 1:
            basis = sc_basis;
            basis_size = 1;
            break;

        case 2:
            basis = bcc_basis;
            basis_size = 2;
            break;

        case 3:
            basis = fcc_basis;
            basis_size = 4;
            break;

        case 4:
            basis = diamond_basis;
            basis_size = 8;
            break;

        default:
            return 0;
    }

    size_t required_count =
        periodic_atom_count(
            structure_choice,
            nx,
            ny,
            nz
        );

    if (capacity < required_count) {
        return 0;
    }

    size_t atom_index = 0;

    // These loops use < rather than <=.
    //
    // This prevents duplicated atoms from appearing on opposite
    // periodic boundaries.
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz; iz++) {
                for (int basis_index = 0;
                     basis_index < basis_size;
                     basis_index++) {

                    snprintf(
                        atoms[atom_index].element,
                        sizeof(atoms[atom_index].element),
                        "%s",
                        element
                    );

                    atoms[atom_index].x =
                        (
                            ix
                            + basis[basis_index][0]
                        )
                        * lattice_constant;

                    atoms[atom_index].y =
                        (
                            iy
                            + basis[basis_index][1]
                        )
                        * lattice_constant;

                    atoms[atom_index].z =
                        (
                            iz
                            + basis[basis_index][2]
                        )
                        * lattice_constant;

                    atom_index++;
                }
            }
        }
    }

    return atom_index;
}

// Calculate a minimum-image displacement vector.
//
// The vector points from atom i to atom j:
//
// displacement = position_j - position_i
void energy_pbc_displacement(
    const Atom *atom_i,
    const Atom *atom_j,
    double lx,
    double ly,
    double lz,
    double displacement[3]
)
{
    displacement[0] =
        atom_j->x - atom_i->x;

    displacement[1] =
        atom_j->y - atom_i->y;

    displacement[2] =
        atom_j->z - atom_i->z;

    displacement[0] -=
        lx * round(displacement[0] / lx);

    displacement[1] -=
        ly * round(displacement[1] / ly);

    displacement[2] -=
        lz * round(displacement[2] / lz);
}

// Check whether the potential cutoff is less than half of
// every supercell length.
int minimum_image_cutoff_is_valid(
    double lx,
    double ly,
    double lz,
    double cutoff
)
{
    if (lx <= 0.0 ||
        ly <= 0.0 ||
        lz <= 0.0 ||
        cutoff <= 0.0) {

        return 0;
    }

    double minimum_length = lx;

    if (ly < minimum_length) {
        minimum_length = ly;
    }

    if (lz < minimum_length) {
        minimum_length = lz;
    }

    return 2.0 * cutoff < minimum_length;
}