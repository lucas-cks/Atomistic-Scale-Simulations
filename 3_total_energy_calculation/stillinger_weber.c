#include "energy.h"

#include <math.h>

#define SW_MINIMUM_DISTANCE_SQUARED 1.0e-24

// Original Stillinger-Weber silicon parameters.
//
// epsilon is in eV.
// sigma is in Angstrom.
// All remaining parameters are dimensionless.
const SWParameters SILICON_SW_PARAMETERS = {
    2.1683,              // epsilon
    2.0951,              // sigma
    1.80,                // cutoff_a
    21.0,                // lambda
    1.20,                // gamma
    -0.333333333333,     // cos_theta_0
    7.049556277,         // A
    0.6022245584,        // B
    4.0,                 // p
    0.0                  // q
};

// Calculate the squared length of a vector.
static double sw_vector_length_squared(
    const double vector[3]
)
{
    return vector[0] * vector[0]
         + vector[1] * vector[1]
         + vector[2] * vector[2];
}

// Calculate the dot product of two vectors.
static double sw_dot_product(
    const double vector_1[3],
    const double vector_2[3]
)
{
    return vector_1[0] * vector_2[0]
         + vector_1[1] * vector_2[1]
         + vector_1[2] * vector_2[2];
}

// Calculate one SW two-body contribution.
static double sw_two_body_pair_energy(
    double distance,
    const SWParameters *parameters
)
{
    double cutoff =
        parameters->cutoff_a
        * parameters->sigma;

    if (distance >= cutoff) {
        return 0.0;
    }

    double sigma_over_r =
        parameters->sigma / distance;

    double radial_part =
        parameters->B
        * pow(sigma_over_r, parameters->p)
        - pow(sigma_over_r, parameters->q);

    double cutoff_part =
        exp(
            parameters->sigma
            / (distance - cutoff)
        );

    return parameters->A
         * parameters->epsilon
         * radial_part
         * cutoff_part;
}

// Calculate one SW three-body contribution.
//
// Atom i is the central atom.
// vector_ij points from atom i to atom j.
// vector_ik points from atom i to atom k.
static double sw_three_body_angle_energy(
    const double vector_ij[3],
    double distance_ij,
    const double vector_ik[3],
    double distance_ik,
    const SWParameters *parameters
)
{
    double cutoff =
        parameters->cutoff_a
        * parameters->sigma;

    if (distance_ij >= cutoff ||
        distance_ik >= cutoff) {

        return 0.0;
    }

    double cosine =
        sw_dot_product(vector_ij, vector_ik)
        / (distance_ij * distance_ik);

    // Protect against small floating-point errors.
    if (cosine > 1.0) {
        cosine = 1.0;
    } else if (cosine < -1.0) {
        cosine = -1.0;
    }

    double angular_difference =
        cosine - parameters->cos_theta_0;

    double cutoff_part_ij =
        exp(
            parameters->gamma
            * parameters->sigma
            / (distance_ij - cutoff)
        );

    double cutoff_part_ik =
        exp(
            parameters->gamma
            * parameters->sigma
            / (distance_ik - cutoff)
        );

    return parameters->lambda
         * parameters->epsilon
         * angular_difference
         * angular_difference
         * cutoff_part_ij
         * cutoff_part_ik;
}

// Calculate the total Stillinger-Weber energy.
double calculate_sw_total_energy(
    const Atom atoms[],
    size_t number_of_atoms,
    double lx,
    double ly,
    double lz,
    const SWParameters *parameters,
    double *two_body_energy,
    double *three_body_energy
)
{
    if (two_body_energy != NULL) {
        *two_body_energy = NAN;
    }

    if (three_body_energy != NULL) {
        *three_body_energy = NAN;
    }

    if (atoms == NULL ||
        parameters == NULL ||
        number_of_atoms == 0 ||
        lx <= 0.0 ||
        ly <= 0.0 ||
        lz <= 0.0 ||
        parameters->epsilon <= 0.0 ||
        parameters->sigma <= 0.0 ||
        parameters->cutoff_a <= 0.0) {

        return NAN;
    }

    double pair_total = 0.0;
    double angle_total = 0.0;

    double cutoff =
        parameters->cutoff_a
        * parameters->sigma;

    double cutoff_squared =
        cutoff * cutoff;

    // Calculate the two-body energy.
    //
    // Every pair is counted once because j begins at i + 1.
    for (size_t i = 0; i < number_of_atoms; i++) {
        for (size_t j = i + 1;
             j < number_of_atoms;
             j++) {

            double vector_ij[3];

            energy_pbc_displacement(
                &atoms[i],
                &atoms[j],
                lx,
                ly,
                lz,
                vector_ij
            );

            double distance_squared =
                sw_vector_length_squared(vector_ij);

            if (distance_squared <
                SW_MINIMUM_DISTANCE_SQUARED) {

                return NAN;
            }

            if (distance_squared < cutoff_squared) {
                double distance =
                    sqrt(distance_squared);

                pair_total +=
                    sw_two_body_pair_energy(
                        distance,
                        parameters
                    );
            }
        }
    }

    // Calculate the three-body energy.
    //
    // Atom i is the central atom of angle j-i-k.
    //
    // For each central atom, k begins at j + 1 so every
    // unordered neighbor pair is counted once.
    for (size_t i = 0; i < number_of_atoms; i++) {
        for (size_t j = 0; j < number_of_atoms; j++) {
            if (j == i) {
                continue;
            }

            double vector_ij[3];

            energy_pbc_displacement(
                &atoms[i],
                &atoms[j],
                lx,
                ly,
                lz,
                vector_ij
            );

            double distance_ij_squared =
                sw_vector_length_squared(vector_ij);

            if (distance_ij_squared <
                SW_MINIMUM_DISTANCE_SQUARED) {

                return NAN;
            }

            if (distance_ij_squared >= cutoff_squared) {
                continue;
            }

            double distance_ij =
                sqrt(distance_ij_squared);

            for (size_t k = j + 1;
                 k < number_of_atoms;
                 k++) {

                if (k == i) {
                    continue;
                }

                double vector_ik[3];

                energy_pbc_displacement(
                    &atoms[i],
                    &atoms[k],
                    lx,
                    ly,
                    lz,
                    vector_ik
                );

                double distance_ik_squared =
                    sw_vector_length_squared(vector_ik);

                if (distance_ik_squared <
                    SW_MINIMUM_DISTANCE_SQUARED) {

                    return NAN;
                }

                if (distance_ik_squared >=
                    cutoff_squared) {

                    continue;
                }

                double distance_ik =
                    sqrt(distance_ik_squared);

                angle_total +=
                    sw_three_body_angle_energy(
                        vector_ij,
                        distance_ij,
                        vector_ik,
                        distance_ik,
                        parameters
                    );
            }
        }
    }

    if (two_body_energy != NULL) {
        *two_body_energy = pair_total;
    }

    if (three_body_energy != NULL) {
        *three_body_energy = angle_total;
    }

    return pair_total + angle_total;
}