#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One atom and its information */
typedef struct {
    char element[3];
    double x;
    double y;
    double z;
} Atom;


/* Function declarations */
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

int write_xyz(
    const char *filename,
    const Atom *atoms,
    int number_of_atoms,
    const char *comment
);


/* Generate a simple-cubic structure */
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

    /* Number of unique lattice points in a finite block */
    *number_of_atoms = (nx + 1) * (ny + 1) * (nz + 1);

    /* Allocate enough memory for all atoms */
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    /* Generate every simple-cubic lattice point */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }
}

/* Generate a body-centred cubic structure */
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

    /* Number of unique lattice points in a finite block */
    *number_of_atoms = (nx + 1) * (ny + 1) * (nz + 1) + nx * ny * nz;

    /* Allocate enough memory for all atoms */
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    /* Generate bcc corner lattice point */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }

    /* Generate body-centre atom inside each unit cell */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = (i + 0.5) * lattice_constant;
                (*atoms)[idx].y = (j + 0.5) * lattice_constant;
                (*atoms)[idx].z = (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }
}

/* Generate a face-centred cubic structure */
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

    /* Number of unique lattice points in a finite block */
    *number_of_atoms =
        (nx + 1) * (ny + 1) * (nz + 1)  /* Corners */
        + nx * ny * (nz + 1)             /* xy faces */
        + nx * (ny + 1) * nz             /* xz faces */
        + (nx + 1) * ny * nz;            /* yz faces */
    /* Allocate enough memory for all atoms */
    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    /* Generate fcc corner lattice point */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }

    /* Face centers parallel to the xy-plane */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = (i + 0.5) * lattice_constant;
                (*atoms)[idx].y = (j + 0.5) * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }

    /* Face centers parallel to the xz-plane */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = (i + 0.5) * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    /* Face centers parallel to the yz-plane */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = (j + 0.5) * lattice_constant;
                (*atoms)[idx].z = (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }
}

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

    /*
     * Diamond cubic:
     *
     * 1. FCC corner and face-center atoms
     * 2. Four additional internal atoms per unit cell
     */
    *number_of_atoms =
        (nx + 1) * (ny + 1) * (nz + 1)  /* Corners */
        + nx * ny * (nz + 1)             /* xy faces */
        + nx * (ny + 1) * nz             /* xz faces */
        + (nx + 1) * ny * nz             /* yz faces */
        + 4 * nx * ny * nz;              /* Internal atoms */

    *atoms = malloc(
        (size_t)(*number_of_atoms) * sizeof **atoms
    );

    if (*atoms == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    /* Generate corner atoms */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }

    /* Generate face centers parallel to the xy-plane */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k <= nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = (i + 0.5) * lattice_constant;
                (*atoms)[idx].y = (j + 0.5) * lattice_constant;
                (*atoms)[idx].z = k * lattice_constant;

                idx++;
            }
        }
    }

    /* Generate face centers parallel to the xz-plane */
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j <= ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = (i + 0.5) * lattice_constant;
                (*atoms)[idx].y = j * lattice_constant;
                (*atoms)[idx].z = (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    /* Generate face centers parallel to the yz-plane */
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                strcpy((*atoms)[idx].element, element);

                (*atoms)[idx].x = i * lattice_constant;
                (*atoms)[idx].y = (j + 0.5) * lattice_constant;
                (*atoms)[idx].z = (k + 0.5) * lattice_constant;

                idx++;
            }
        }
    }

    /*
     * Four internal atoms belonging to the second,
     * displaced FCC lattice.
     */
    const double internal_basis[4][3] = {
        {0.25, 0.25, 0.25},
        {0.25, 0.75, 0.75},
        {0.75, 0.25, 0.75},
        {0.75, 0.75, 0.25}
    };

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {

                for (int b = 0; b < 4; b++) {

                    strcpy((*atoms)[idx].element, element);

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


// Main function
int main(void)
{
    char element[3];
    char title[96];
    char filename[100];
    char comment[200];
    int structure_choice;

    double lattice_constant;
    int nx;
    int ny;
    int nz;

    Atom *atoms = NULL;
    int number_of_atoms = 0;

    printf("Cubic Crystal Generator\n");


    printf ("Choose the crystal structure:\n");
    printf ("1. Simple Cubic\n");
    printf ("2. Body-Centred Cubic\n");
    printf ("3. Face-Centred Cubic\n");
    printf ("4. Diamond Cubic\n");
    printf("Enter your choice (1-4): ");

    if (scanf("%d", &structure_choice) != 1 ||
    structure_choice < 1 ||
    structure_choice > 4) {
    printf("Invalid structure choice.\n");
    return 1;
}

    /* Read the chemical element */
    printf("Enter the element symbol (e.g.: Si): ");

    if (scanf("%2s", element) != 1) {
        printf("Invalid element.\n");
        return 1;
    }

    /* Read the lattice constant */
    printf("Enter the lattice constant in Angstroms: ");

    if (scanf("%lf", &lattice_constant) != 1 ||
        lattice_constant <= 0.0) {
        printf("Invalid lattice constant.\n");
        return 1;
    }

    /* Read the number of unit cells */
    printf("Enter the number of unit cells in the x-direction: ");

    if (scanf("%d", &nx) != 1 || nx < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    printf("Enter the number of unit cells in the y-direction: ");

    if (scanf("%d", &ny) != 1 || ny < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    printf("Enter the number of unit cells in the z-direction: ");

    if (scanf("%d", &nz) != 1 || nz < 1) {
        printf("Invalid number of unit cells.\n");
        return 1;
    }

    /* Read the title*/
    printf("Enter the file title (without .xyz): ");

    if (scanf("%95s", title) != 1) {
        printf("Invalid filename.\n");
        return 1;
    }

    /* Read a comment; spaces are permitted */
    printf("Enter a comment for the XYZ file: ");

    if (scanf(" %199[^\n]", comment) != 1) {
        printf("Invalid comment.\n");
        return 1;
    }

    /* Allocate memory and generate the atom positions */

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

            snprintf(filename, sizeof(filename), "%s_sc.xyz", title);
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

            snprintf(filename, sizeof(filename), "%s_bcc.xyz", title);
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

            snprintf(filename, sizeof(filename), "%s_fcc.xyz", title);
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

            snprintf(filename, sizeof(filename), "%s_diamond.xyz", title);
            break;

        default:
            printf("Invalid structure choice.\n");
            return 1;
    }

    /* Write the structure to the XYZ file */
    if (write_xyz(
            filename,
            atoms,
            number_of_atoms,
            comment
        ) != 0) {

        free(atoms);
        return 1;
    }

    printf("\nSuccessfully generated \"%s\".\n", filename);
    printf("Number of atoms: %d\n", number_of_atoms);

    /* Release the memory created by malloc() */
    free(atoms);
    atoms = NULL;

    return 0;
}

/*
  Write the atom array to an XYZ file.
 
  XYZ format:
 
  number of atoms
  comment
  element x y z
  element x y z
  ...
*/
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

    /* First line: number of atoms */
    fprintf(output, "%d\n", number_of_atoms);

    /* Second line: comment */
    fprintf(output, "%s\n", comment);

    /* Remaining lines: element and coordinates */
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