#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_NEIGHBORS 32

// One atom and its information
typedef struct {
    char element[3];
    double x;
    double y;
    double z;
} Atom;

// One neighbor and its information
typedef struct {
    int atom_index;
    int image[3];
    double vector[3];
    double distance;
} Neighbor;

// Neighbor list belonging to one atom
typedef struct {
    int count;
    Neighbor neighbors[MAX_NEIGHBORS];
} NeighborList;

// Primitive vectors for the crystal structure
double primitive_vectors[3][3] = {0.0};

// Reciprocal vectors for the crystal structure
double reciprocal_vectors[3][3] = {0.0};


// Function declarations

void generate_sc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
);

void generate_bcc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
);

void generate_fcc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
);

void generate_diamond(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
);

void set_primitive_vectors_sc(
    double primitive_vectors[3][3],
    double lattice_constant
);

void set_primitive_vectors_bcc(
    double primitive_vectors[3][3],
    double lattice_constant
);

void set_primitive_vectors_fcc_diamond(
    double primitive_vectors[3][3],
    double lattice_constant
);

void print_primitive_vectors(
    const double primitive_vectors[3][3]
);

double primitive_cell_volume(
    const double primitive_vectors[3][3]
);

void set_reciprocal_vectors(
    double reciprocal_vectors[3][3]
);

void print_reciprocal_vectors(
    const double reciprocal_vectors[3][3]
);

void vector_under_pbc(
    double x1,
    double y1,
    double z1,
    double fractional_pbc[3],
    double cartesian_pbc[3]
);

void vector_from_two_atoms(
    const Atom *atoms,
    int atom1,
    int atom2,
    double original_cartesian[3],
    double fractional_pbc[3],
    double cartesian_pbc[3]
);

int set_primitive_cell_atoms(
    Atom unit_cell_atoms[2],
    int structure_choice,
    const char *element,
    double lattice_constant
);

double first_neighbor_distance(
    int structure_choice,
    double lattice_constant
);

double second_neighbor_distance(
    int structure_choice,
    double lattice_constant
);

int build_neighbor_lists(
    const Atom *unit_cell_atoms,
    int number_of_unit_cell_atoms,
    double cutoff,
    NeighborList *neighbor_lists
);

int write_neighbor_lists(
    const char *filename,
    const Atom *unit_cell_atoms,
    int number_of_unit_cell_atoms,
    const NeighborList *neighbor_lists,
    double cutoff
);

int write_xyz(
    const char *filename,
    const Atom *atoms,
    int number_of_atoms,
    const char *comment
);


// Generate a simple-cubic structure

void generate_sc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
)
{
    int idx = 0;

    // Number of unique lattice points in a finite block
    *number_of_atoms =
        (nx + 1) * (ny + 1) * (nz + 1);

    // Allocate enough memory for all atoms
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // Generate every simple-cubic lattice point
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }
}


// Generate a body-centred cubic structure

void generate_bcc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
)
{
    int idx = 0;

    // Number of unique lattice points in a finite block
    *number_of_atoms =
          (nx + 1) * (ny + 1) * (nz + 1)
        + nx * ny * nz;

    // Allocate enough memory for all atoms
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // Generate BCC corner lattice points
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }

    // Generate the body-centre atom inside each unit cell
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    (i + 0.5) * lattice_constant;

                (*atoms)[idx].y =
                    (j + 0.5) * lattice_constant;

                (*atoms)[idx].z =
                    (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }
}


// Generate a face-centred cubic structure

void generate_fcc(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
)
{
    int idx = 0;

    // Number of unique lattice points in a finite block
    *number_of_atoms =
          (nx + 1) * (ny + 1) * (nz + 1)
        + nx * ny * (nz + 1)
        + nx * (ny + 1) * nz
        + (nx + 1) * ny * nz;

    // Allocate enough memory for all atoms
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // Generate FCC corner lattice points
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the xy-plane
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    (i + 0.5) * lattice_constant;

                (*atoms)[idx].y =
                    (j + 0.5) * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the xz-plane
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    (i + 0.5) * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the yz-plane
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    (j + 0.5) * lattice_constant;

                (*atoms)[idx].z =
                    (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }
}


// Generate a diamond-cubic structure

void generate_diamond(
    Atom **atoms,
    int *number_of_atoms,
    const char *element,
    double lattice_constant,
    int nx,
    int ny,
    int nz
)
{
    int idx = 0;

    // Diamond contains FCC lattice points and four internal atoms
    *number_of_atoms =
          (nx + 1) * (ny + 1) * (nz + 1)
        + nx * ny * (nz + 1)
        + nx * (ny + 1) * nz
        + (nx + 1) * ny * nz
        + 4 * nx * ny * nz;

    // Allocate enough memory for all atoms
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // Generate corner atoms
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the xy-plane
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    (i + 0.5) * lattice_constant;

                (*atoms)[idx].y =
                    (j + 0.5) * lattice_constant;

                (*atoms)[idx].z =
                    k * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the xz-plane
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    (i + 0.5) * lattice_constant;

                (*atoms)[idx].y =
                    j * lattice_constant;

                (*atoms)[idx].z =
                    (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    // Generate face centres parallel to the yz-plane
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x =
                    i * lattice_constant;

                (*atoms)[idx].y =
                    (j + 0.5) * lattice_constant;

                (*atoms)[idx].z =
                    (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    // Fractional positions of the four internal atoms
    const double internal_basis[4][3] = {
        {0.25, 0.25, 0.25},
        {0.25, 0.75, 0.75},
        {0.75, 0.25, 0.75},
        {0.75, 0.75, 0.25}
    };

    // Generate the internal atoms
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                for (int b = 0; b < 4; b++) {

                    strcpy(
                        (*atoms)[idx].element,
                        element
                    );

                    (*atoms)[idx].x =
                        (i + internal_basis[b][0])
                        * lattice_constant;

                    (*atoms)[idx].y =
                        (j + internal_basis[b][1])
                        * lattice_constant;

                    (*atoms)[idx].z =
                        (k + internal_basis[b][2])
                        * lattice_constant;

                    idx++;
                }
            }
        }
    }
}


// Set the simple-cubic primitive vectors

void set_primitive_vectors_sc(
    double primitive_vectors[3][3],
    double lattice_constant
)
{
    primitive_vectors[0][0] = lattice_constant;
    primitive_vectors[0][1] = 0.0;
    primitive_vectors[0][2] = 0.0;

    primitive_vectors[1][0] = 0.0;
    primitive_vectors[1][1] = lattice_constant;
    primitive_vectors[1][2] = 0.0;

    primitive_vectors[2][0] = 0.0;
    primitive_vectors[2][1] = 0.0;
    primitive_vectors[2][2] = lattice_constant;
}


// Set the body-centred cubic primitive vectors

void set_primitive_vectors_bcc(
    double primitive_vectors[3][3],
    double lattice_constant
)
{
    primitive_vectors[0][0] =
        -lattice_constant / 2.0;

    primitive_vectors[0][1] =
         lattice_constant / 2.0;

    primitive_vectors[0][2] =
         lattice_constant / 2.0;

    primitive_vectors[1][0] =
         lattice_constant / 2.0;

    primitive_vectors[1][1] =
        -lattice_constant / 2.0;

    primitive_vectors[1][2] =
         lattice_constant / 2.0;

    primitive_vectors[2][0] =
         lattice_constant / 2.0;

    primitive_vectors[2][1] =
         lattice_constant / 2.0;

    primitive_vectors[2][2] =
        -lattice_constant / 2.0;
}


// Set the FCC and diamond primitive vectors

void set_primitive_vectors_fcc_diamond(
    double primitive_vectors[3][3],
    double lattice_constant
)
{
    primitive_vectors[0][0] = 0.0;

    primitive_vectors[0][1] =
        lattice_constant / 2.0;

    primitive_vectors[0][2] =
        lattice_constant / 2.0;

    primitive_vectors[1][0] =
        lattice_constant / 2.0;

    primitive_vectors[1][1] = 0.0;

    primitive_vectors[1][2] =
        lattice_constant / 2.0;

    primitive_vectors[2][0] =
        lattice_constant / 2.0;

    primitive_vectors[2][1] =
        lattice_constant / 2.0;

    primitive_vectors[2][2] = 0.0;
}


// Print the primitive vectors

void print_primitive_vectors(
    const double primitive_vectors[3][3]
)
{
    printf("\nPrimitive lattice vectors:\n");

    for (int i = 0; i < 3; i++) {
        printf(
            "a%d = (%.5f, %.5f, %.5f) Angstrom\n",
            i + 1,
            primitive_vectors[i][0],
            primitive_vectors[i][1],
            primitive_vectors[i][2]
        );
    }
}


// Calculate the primitive-cell volume

double primitive_cell_volume(
    const double primitive_vectors[3][3]
)
{
    double determinant =
          primitive_vectors[0][0]
        * (
              primitive_vectors[1][1]
            * primitive_vectors[2][2]
            - primitive_vectors[1][2]
            * primitive_vectors[2][1]
          )

        - primitive_vectors[0][1]
        * (
              primitive_vectors[1][0]
            * primitive_vectors[2][2]
            - primitive_vectors[1][2]
            * primitive_vectors[2][0]
          )

        + primitive_vectors[0][2]
        * (
              primitive_vectors[1][0]
            * primitive_vectors[2][1]
            - primitive_vectors[1][1]
            * primitive_vectors[2][0]
          );

    return fabs(determinant);
}


// Calculate the reciprocal lattice vectors

void set_reciprocal_vectors(
    double reciprocal_vectors[3][3]
)
{
    double volume =
        primitive_cell_volume(primitive_vectors);

    // First reciprocal vector
    reciprocal_vectors[0][0] =
        (
              primitive_vectors[1][1]
            * primitive_vectors[2][2]
            - primitive_vectors[1][2]
            * primitive_vectors[2][1]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[0][1] =
        (
              primitive_vectors[1][2]
            * primitive_vectors[2][0]
            - primitive_vectors[1][0]
            * primitive_vectors[2][2]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[0][2] =
        (
              primitive_vectors[1][0]
            * primitive_vectors[2][1]
            - primitive_vectors[1][1]
            * primitive_vectors[2][0]
        ) * 2.0 * M_PI / volume;

    // Second reciprocal vector
    reciprocal_vectors[1][0] =
        (
             -primitive_vectors[0][1]
            * primitive_vectors[2][2]
            + primitive_vectors[0][2]
            * primitive_vectors[2][1]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[1][1] =
        (
             -primitive_vectors[0][2]
            * primitive_vectors[2][0]
            + primitive_vectors[0][0]
            * primitive_vectors[2][2]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[1][2] =
        (
             -primitive_vectors[0][0]
            * primitive_vectors[2][1]
            + primitive_vectors[0][1]
            * primitive_vectors[2][0]
        ) * 2.0 * M_PI / volume;

    // Third reciprocal vector
    reciprocal_vectors[2][0] =
        (
              primitive_vectors[0][1]
            * primitive_vectors[1][2]
            - primitive_vectors[0][2]
            * primitive_vectors[1][1]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[2][1] =
        (
              primitive_vectors[0][2]
            * primitive_vectors[1][0]
            - primitive_vectors[0][0]
            * primitive_vectors[1][2]
        ) * 2.0 * M_PI / volume;

    reciprocal_vectors[2][2] =
        (
              primitive_vectors[0][0]
            * primitive_vectors[1][1]
            - primitive_vectors[0][1]
            * primitive_vectors[1][0]
        ) * 2.0 * M_PI / volume;
}


// Print the reciprocal lattice vectors

void print_reciprocal_vectors(
    const double reciprocal_vectors[3][3]
)
{
    printf("\nReciprocal lattice vectors:\n");

    for (int i = 0; i < 3; i++) {
        printf(
            "b%d = (%.5f, %.5f, %.5f) Angstrom^-1\n",
            i + 1,
            reciprocal_vectors[i][0],
            reciprocal_vectors[i][1],
            reciprocal_vectors[i][2]
        );
    }
}


// Apply periodic boundary conditions to a Cartesian vector

void vector_under_pbc(
    double x1,
    double y1,
    double z1,
    double fractional_pbc[3],
    double cartesian_pbc[3]
)
{
    // Convert the Cartesian vector to fractional components
    fractional_pbc[0] =
        (
              x1 * reciprocal_vectors[0][0]
            + y1 * reciprocal_vectors[0][1]
            + z1 * reciprocal_vectors[0][2]
        ) / (2.0 * M_PI);

    fractional_pbc[1] =
        (
              x1 * reciprocal_vectors[1][0]
            + y1 * reciprocal_vectors[1][1]
            + z1 * reciprocal_vectors[1][2]
        ) / (2.0 * M_PI);

    fractional_pbc[2] =
        (
              x1 * reciprocal_vectors[2][0]
            + y1 * reciprocal_vectors[2][1]
            + z1 * reciprocal_vectors[2][2]
        ) / (2.0 * M_PI);

    // Wrap each component into the interval (-0.5, 0.5]
    for (int i = 0; i < 3; i++) {
        fractional_pbc[i] =
              fractional_pbc[i]
            - ceil(fractional_pbc[i] - 0.5);
    }

    // Convert the wrapped vector back to Cartesian components
    for (int component = 0; component < 3; component++) {
        cartesian_pbc[component] =
              fractional_pbc[0]
            * primitive_vectors[0][component]

            + fractional_pbc[1]
            * primitive_vectors[1][component]

            + fractional_pbc[2]
            * primitive_vectors[2][component];
    }
}


// Calculate the vector between two atoms

void vector_from_two_atoms(
    const Atom *atoms,
    int atom1,
    int atom2,
    double original_cartesian[3],
    double fractional_pbc[3],
    double cartesian_pbc[3]
)
{
    // Calculate the original Cartesian displacement
    original_cartesian[0] =
        atoms[atom2].x - atoms[atom1].x;

    original_cartesian[1] =
        atoms[atom2].y - atoms[atom1].y;

    original_cartesian[2] =
        atoms[atom2].z - atoms[atom1].z;

    // Apply periodic boundary conditions
    vector_under_pbc(
        original_cartesian[0],
        original_cartesian[1],
        original_cartesian[2],
        fractional_pbc,
        cartesian_pbc
    );
}


// Set the atoms belonging to one primitive cell

int set_primitive_cell_atoms(
    Atom unit_cell_atoms[2],
    int structure_choice,
    const char *element,
    double lattice_constant
)
{
    // SC, BCC and FCC have one atom per primitive cell
    strcpy(unit_cell_atoms[0].element, element);

    unit_cell_atoms[0].x = 0.0;
    unit_cell_atoms[0].y = 0.0;
    unit_cell_atoms[0].z = 0.0;

    // Diamond has two atoms per primitive cell
    if (structure_choice == 4) {
        strcpy(unit_cell_atoms[1].element, element);

        unit_cell_atoms[1].x =
            lattice_constant / 4.0;

        unit_cell_atoms[1].y =
            lattice_constant / 4.0;

        unit_cell_atoms[1].z =
            lattice_constant / 4.0;

        return 2;
    }

    return 1;
}


// Return the first-neighbor distance

double first_neighbor_distance(
    int structure_choice,
    double lattice_constant
)
{
    switch (structure_choice) {
        case 1:
            return lattice_constant;

        case 2:
            return
                sqrt(3.0) * lattice_constant / 2.0;

        case 3:
            return
                lattice_constant / sqrt(2.0);

        case 4:
            return
                sqrt(3.0) * lattice_constant / 4.0;

        default:
            return 0.0;
    }
}


// Return the second-neighbor distance

double second_neighbor_distance(
    int structure_choice,
    double lattice_constant
)
{
    switch (structure_choice) {
        case 1:
            return
                sqrt(2.0) * lattice_constant;

        case 2:
            return lattice_constant;

        case 3:
            return lattice_constant;

        case 4:
            return
                lattice_constant / sqrt(2.0);

        default:
            return 0.0;
    }
}


// Build the neighbor list for every primitive-cell atom

int build_neighbor_lists(
    const Atom *unit_cell_atoms,
    int number_of_unit_cell_atoms,
    double cutoff,
    NeighborList *neighbor_lists
)
{
    const double zero_tolerance = 1.0e-8;
    const double cutoff_tolerance = 1.0e-8;

    // Set the initial number of neighbors to zero
    for (int i = 0;
         i < number_of_unit_cell_atoms;
         i++) {

        neighbor_lists[i].count = 0;
    }

    // Examine each atom as the central atom
    for (int i = 0;
         i < number_of_unit_cell_atoms;
         i++) {

        // Examine each basis atom as a possible neighbor
        for (int j = 0;
             j < number_of_unit_cell_atoms;
             j++) {

            // Examine the central cell and surrounding cells
            for (int n1 = -1; n1 <= 1; n1++) {
                for (int n2 = -1; n2 <= 1; n2++) {
                    for (int n3 = -1; n3 <= 1; n3++) {

                        // Exclude the central atom from its own list
                        if (i == j &&
                            n1 == 0 &&
                            n2 == 0 &&
                            n3 == 0) {

                            continue;
                        }

                        // Calculate the Cartesian neighbor vector
                        double dx =
                              unit_cell_atoms[j].x
                            - unit_cell_atoms[i].x
                            + n1 * primitive_vectors[0][0]
                            + n2 * primitive_vectors[1][0]
                            + n3 * primitive_vectors[2][0];

                        double dy =
                              unit_cell_atoms[j].y
                            - unit_cell_atoms[i].y
                            + n1 * primitive_vectors[0][1]
                            + n2 * primitive_vectors[1][1]
                            + n3 * primitive_vectors[2][1];

                        double dz =
                              unit_cell_atoms[j].z
                            - unit_cell_atoms[i].z
                            + n1 * primitive_vectors[0][2]
                            + n2 * primitive_vectors[1][2]
                            + n3 * primitive_vectors[2][2];

                        // Calculate the neighbor distance
                        double distance = sqrt(
                              dx * dx
                            + dy * dy
                            + dz * dz
                        );

                        // Store atoms inside the selected cutoff
                        if (distance > zero_tolerance &&
                            distance <=
                                cutoff + cutoff_tolerance) {

                            int index =
                                neighbor_lists[i].count;

                            if (index >= MAX_NEIGHBORS) {
                                printf(
                                    "Too many neighbors for "
                                    "primitive-cell atom %d.\n",
                                    i + 1
                                );

                                return 0;
                            }

                            Neighbor *neighbor =
                                &neighbor_lists[i]
                                    .neighbors[index];

                            neighbor->atom_index = j;

                            neighbor->image[0] = n1;
                            neighbor->image[1] = n2;
                            neighbor->image[2] = n3;

                            neighbor->vector[0] = dx;
                            neighbor->vector[1] = dy;
                            neighbor->vector[2] = dz;

                            neighbor->distance = distance;

                            neighbor_lists[i].count++;
                        }
                    }
                }
            }
        }
    }

    return 1;
}


// Write the neighbor lists to a text file

int write_neighbor_lists(
    const char *filename,
    const Atom *unit_cell_atoms,
    int number_of_unit_cell_atoms,
    const NeighborList *neighbor_lists,
    double cutoff
)
{
    FILE *output = fopen(filename, "w");

    if (output == NULL) {
        perror("Could not open the neighbor-list file");
        return 1;
    }

    fprintf(
        output,
        "First-nearest-neighbor lists\n"
    );

    fprintf(
        output,
        "Distance cutoff: %.8f Angstrom\n\n",
        cutoff
    );

    // Write the list belonging to each atom
    for (int i = 0;
         i < number_of_unit_cell_atoms;
         i++) {

        fprintf(
            output,
            "Atom %d: %s "
            "(%.8f, %.8f, %.8f)\n",
            i + 1,
            unit_cell_atoms[i].element,
            unit_cell_atoms[i].x,
            unit_cell_atoms[i].y,
            unit_cell_atoms[i].z
        );

        fprintf(
            output,
            "Number of neighbors: %d\n",
            neighbor_lists[i].count
        );

        // Write each individual neighbor
        for (int k = 0;
             k < neighbor_lists[i].count;
             k++) {

            const Neighbor *neighbor =
                &neighbor_lists[i].neighbors[k];

            fprintf(
                output,
                "  Neighbor %d: "
                "atom %d, "
                "image (%+d, %+d, %+d), "
                "vector (%.8f, %.8f, %.8f), "
                "distance %.8f Angstrom\n",
                k + 1,
                neighbor->atom_index + 1,
                neighbor->image[0],
                neighbor->image[1],
                neighbor->image[2],
                neighbor->vector[0],
                neighbor->vector[1],
                neighbor->vector[2],
                neighbor->distance
            );
        }

        fprintf(output, "\n");
    }

    fclose(output);

    return 0;
}


// Write the atom array to an XYZ file

int write_xyz(
    const char *filename,
    const Atom *atoms,
    int number_of_atoms,
    const char *comment
)
{
    FILE *output = fopen(filename, "w");

    if (output == NULL) {
        perror("Could not open the output file");
        return 1;
    }

    // First line contains the number of atoms
    fprintf(output, "%d\n", number_of_atoms);

    // Second line contains the comment
    fprintf(output, "%s\n", comment);

    // Remaining lines contain elements and coordinates
    for (int i = 0; i < number_of_atoms; i++) {
        fprintf(
            output,
            "%s %.5f %.5f %.5f\n",
            atoms[i].element,
            atoms[i].x,
            atoms[i].y,
            atoms[i].z
        );
    }

    fclose(output);

    return 0;
}


// Main function

int main(void)
{
    char element[3];
    char title[96];
    char filename[128];
    char neighbor_filename[128];
    char comment[200];

    int structure_choice;

    double lattice_constant;

    int nx;
    int ny;
    int nz;

    Atom *atoms = NULL;
    int number_of_atoms = 0;

    printf("Cubic Crystal Generator\n");

    printf("Choose the crystal structure:\n");
    printf("1. Simple Cubic\n");
    printf("2. Body-Centred Cubic\n");
    printf("3. Face-Centred Cubic\n");
    printf("4. Diamond Cubic\n");
    printf("Enter your choice (1-4): ");

    if (scanf("%d", &structure_choice) != 1 ||
        structure_choice < 1 ||
        structure_choice > 4) {

        printf("Invalid structure choice.\n");
        return 1;
    }

    // Read the chemical element
    printf("Enter the element symbol (e.g.: Si): ");

    if (scanf("%2s", element) != 1) {
        printf("Invalid element.\n");
        return 1;
    }

    // Read the lattice constant
    printf(
        "Enter the lattice constant in Angstroms: "
    );

    if (scanf("%lf", &lattice_constant) != 1 ||
        lattice_constant <= 0.0) {

        printf("Invalid lattice constant.\n");
        return 1;
    }

    // Read the number of unit cells
    printf(
        "Enter the number of unit cells "
        "in the x-direction: "
    );

    if (scanf("%d", &nx) != 1 || nx < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    printf(
        "Enter the number of unit cells "
        "in the y-direction: "
    );

    if (scanf("%d", &ny) != 1 || ny < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    printf(
        "Enter the number of unit cells "
        "in the z-direction: "
    );

    if (scanf("%d", &nz) != 1 || nz < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    // Read the output title
    printf(
        "Enter the file title without an extension: "
    );

    if (scanf("%95s", title) != 1) {
        printf("Invalid filename.\n");
        return 1;
    }

    // Read the XYZ comment
    printf("Enter a comment for the XYZ file: ");

    if (scanf(" %199[^\n]", comment) != 1) {
        printf("Invalid comment.\n");
        return 1;
    }

    // Generate the selected crystal structure
    switch (structure_choice) {
        case 1:
            generate_sc(
                &atoms,
                &number_of_atoms,
                element,
                lattice_constant,
                nx,
                ny,
                nz
            );

            set_primitive_vectors_sc(
                primitive_vectors,
                lattice_constant
            );

            snprintf(
                filename,
                sizeof filename,
                "%s_sc.xyz",
                title
            );

            break;

        case 2:
            generate_bcc(
                &atoms,
                &number_of_atoms,
                element,
                lattice_constant,
                nx,
                ny,
                nz
            );

            set_primitive_vectors_bcc(
                primitive_vectors,
                lattice_constant
            );

            snprintf(
                filename,
                sizeof filename,
                "%s_bcc.xyz",
                title
            );

            break;

        case 3:
            generate_fcc(
                &atoms,
                &number_of_atoms,
                element,
                lattice_constant,
                nx,
                ny,
                nz
            );

            set_primitive_vectors_fcc_diamond(
                primitive_vectors,
                lattice_constant
            );

            snprintf(
                filename,
                sizeof filename,
                "%s_fcc.xyz",
                title
            );

            break;

        case 4:
            generate_diamond(
                &atoms,
                &number_of_atoms,
                element,
                lattice_constant,
                nx,
                ny,
                nz
            );

            set_primitive_vectors_fcc_diamond(
                primitive_vectors,
                lattice_constant
            );

            snprintf(
                filename,
                sizeof filename,
                "%s_diamond.xyz",
                title
            );

            break;

        default:
            printf("Invalid structure choice.\n");
            return 1;
    }

    // Calculate and print the lattice vectors
    set_reciprocal_vectors(reciprocal_vectors);

    print_primitive_vectors(primitive_vectors);

    double primitive_volume =
        primitive_cell_volume(primitive_vectors);

    printf(
        "Primitive-cell volume: "
        "%.5f Angstrom^3\n",
        primitive_volume
    );

    print_reciprocal_vectors(reciprocal_vectors);

    // Read two atom numbers for the PBC vector calculation
    int atom_number_1;
    int atom_number_2;

    printf(
        "\nChoose the first atom number (1-%d): ",
        number_of_atoms
    );

    if (scanf("%d", &atom_number_1) != 1 ||
        atom_number_1 < 1 ||
        atom_number_1 > number_of_atoms) {

        printf("Invalid atom number.\n");
        free(atoms);
        return 1;
    }

    printf(
        "Choose the second atom number (1-%d): ",
        number_of_atoms
    );

    if (scanf("%d", &atom_number_2) != 1 ||
        atom_number_2 < 1 ||
        atom_number_2 > number_of_atoms) {

        printf("Invalid atom number.\n");
        free(atoms);
        return 1;
    }

    // Convert atom numbers to zero-based indexes
    int atom1 = atom_number_1 - 1;
    int atom2 = atom_number_2 - 1;

    double original_cartesian[3];
    double fractional_pbc[3];
    double cartesian_pbc[3];

    // Calculate the vector between the selected atoms
    vector_from_two_atoms(
        atoms,
        atom1,
        atom2,
        original_cartesian,
        fractional_pbc,
        cartesian_pbc
    );

    // Print the selected atom positions
    printf(
        "\nAtom %d position: "
        "(%.5f, %.5f, %.5f) Angstrom\n",
        atom_number_1,
        atoms[atom1].x,
        atoms[atom1].y,
        atoms[atom1].z
    );

    printf(
        "Atom %d position: "
        "(%.5f, %.5f, %.5f) Angstrom\n",
        atom_number_2,
        atoms[atom2].x,
        atoms[atom2].y,
        atoms[atom2].z
    );

    // Print the original Cartesian vector
    printf(
        "\nOriginal Cartesian vector "
        "from atom %d to atom %d:\n",
        atom_number_1,
        atom_number_2
    );

    printf(
        "(%.5f, %.5f, %.5f) Angstrom\n",
        original_cartesian[0],
        original_cartesian[1],
        original_cartesian[2]
    );

    // Print the fractional vector under PBC
    printf(
        "\nFractional vector under PBC:\n"
        "(%.5f, %.5f, %.5f)\n",
        fractional_pbc[0],
        fractional_pbc[1],
        fractional_pbc[2]
    );

    // Print the Cartesian vector under PBC
    printf(
        "\nCartesian vector under PBC:\n"
        "(%.5f, %.5f, %.5f) Angstrom\n",
        cartesian_pbc[0],
        cartesian_pbc[1],
        cartesian_pbc[2]
    );

    // Calculate the original distance
    double original_distance = sqrt(
          original_cartesian[0]
        * original_cartesian[0]

        + original_cartesian[1]
        * original_cartesian[1]

        + original_cartesian[2]
        * original_cartesian[2]
    );

    // Calculate the distance under PBC
    double pbc_distance = sqrt(
          cartesian_pbc[0]
        * cartesian_pbc[0]

        + cartesian_pbc[1]
        * cartesian_pbc[1]

        + cartesian_pbc[2]
        * cartesian_pbc[2]
    );

    printf(
        "\nOriginal distance: %.5f Angstrom\n",
        original_distance
    );

    printf(
        "Distance under PBC: %.5f Angstrom\n",
        pbc_distance
    );

    // Create the atoms belonging to one primitive cell
    Atom unit_cell_atoms[2];

    int number_of_unit_cell_atoms =
        set_primitive_cell_atoms(
            unit_cell_atoms,
            structure_choice,
            element,
            lattice_constant
        );

    // Calculate the first- and second-neighbor distances
    double expected_first_distance =
        first_neighbor_distance(
            structure_choice,
            lattice_constant
        );

    double expected_second_distance =
        second_neighbor_distance(
            structure_choice,
            lattice_constant
        );

    double recommended_cutoff =
        1.05 * expected_first_distance;

    printf(
        "\nExpected first-neighbor distance: "
        "%.5f Angstrom\n",
        expected_first_distance
    );

    printf(
        "Expected second-neighbor distance: "
        "%.5f Angstrom\n",
        expected_second_distance
    );

    printf(
        "Recommended first-neighbor cutoff: "
        "%.5f Angstrom\n",
        recommended_cutoff
    );

    // Read the neighbor cutoff
    double neighbor_cutoff;

    printf(
        "Enter the neighbor distance cutoff "
        "in Angstroms: "
    );

    if (scanf("%lf", &neighbor_cutoff) != 1 ||
        neighbor_cutoff <= 0.0) {

        printf("Invalid neighbor cutoff.\n");
        free(atoms);
        return 1;
    }

    // Check whether the cutoff includes first neighbors
    if (neighbor_cutoff <
        expected_first_distance - 1.0e-8) {

        printf(
            "Warning: the cutoff is below the "
            "first-neighbor distance.\n"
        );

        printf(
            "The neighbor lists may be empty.\n"
        );
    }

    // Check whether the cutoff includes second neighbors
    if (neighbor_cutoff >=
        expected_second_distance - 1.0e-8) {

        printf(
            "Warning: the cutoff reaches the "
            "second-neighbor distance.\n"
        );

        printf(
            "The lists may contain more than "
            "first-nearest neighbors.\n"
        );
    }

    // Allocate the neighbor-list arrays
    NeighborList *neighbor_lists = calloc(
        (size_t)number_of_unit_cell_atoms,
        sizeof *neighbor_lists
    );

    if (neighbor_lists == NULL) {
        printf(
            "Memory allocation for neighbor lists "
            "failed.\n"
        );

        free(atoms);
        return 1;
    }

    // Build the neighbor lists
    if (!build_neighbor_lists(
            unit_cell_atoms,
            number_of_unit_cell_atoms,
            neighbor_cutoff,
            neighbor_lists
        )) {

        free(neighbor_lists);
        free(atoms);
        return 1;
    }

    // Create the neighbor-list filename
    snprintf(
        neighbor_filename,
        sizeof neighbor_filename,
        "%s_neighbors.txt",
        title
    );

    // Write the neighbor lists
    if (write_neighbor_lists(
            neighbor_filename,
            unit_cell_atoms,
            number_of_unit_cell_atoms,
            neighbor_lists,
            neighbor_cutoff
        ) != 0) {

        free(neighbor_lists);
        free(atoms);
        return 1;
    }

    printf(
        "\nSuccessfully generated \"%s\".\n",
        neighbor_filename
    );

    // Print the number of neighbors found for each atom
    for (int i = 0;
         i < number_of_unit_cell_atoms;
         i++) {

        printf(
            "Primitive-cell atom %d has "
            "%d neighbors.\n",
            i + 1,
            neighbor_lists[i].count
        );
    }

    // Release the neighbor-list memory
    free(neighbor_lists);
    neighbor_lists = NULL;

    // Write the complete structure to the XYZ file
    if (write_xyz(
            filename,
            atoms,
            number_of_atoms,
            comment
        ) != 0) {

        free(atoms);
        return 1;
    }

    printf(
        "\nSuccessfully generated \"%s\".\n",
        filename
    );

    printf(
        "Number of atoms: %d\n",
        number_of_atoms
    );

    // Release the atom-array memory
    free(atoms);
    atoms = NULL;

    return 0;
}