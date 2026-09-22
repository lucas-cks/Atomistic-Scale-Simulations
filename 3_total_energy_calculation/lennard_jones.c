#include "energy.h"

#include <math.h>

#define LJ_MINIMUM_DISTANCE_SQUARED 1.0e-24

double calculate_lj_total_energy(
    const Atom atoms[],
    size_t number_of_atoms,
    double lx,
    double ly,
    double lz,
    double epsilon,
    double sigma,
    double cutoff
)
{
    if (atoms == NULL ||
        number_of_atoms == 0 ||
        lx <= 0.0 ||
        ly <= 0.0 ||
        lz <= 0.0 ||
        epsilon <= 0.0 ||
        sigma <= 0.0 ||
        cutoff <= 0.0) {

        return NAN;
    }

    double total_energy = 0.0;
    double cutoff_squared = cutoff * cutoff;
    double sigma_squared = sigma * sigma;

    // Each atom pair is counted once because j begins at i + 1.
    for (size_t i = 0; i < number_of_atoms; i++) {
        for (size_t j = i + 1;
             j < number_of_atoms;
             j++) {

            double displacement[3];

            energy_pbc_displacement(
                &atoms[i],
                &atoms[j],
                lx,
                ly,
                lz,
                displacement
            );

            double distance_squared =
                displacement[0] * displacement[0]
                + displacement[1] * displacement[1]
                + displacement[2] * displacement[2];

            // A zero PBC distance normally means that duplicate
            // periodic-boundary atoms were supplied.
            if (distance_squared <
                LJ_MINIMUM_DISTANCE_SQUARED) {

                return NAN;
            }

            if (distance_squared < cutoff_squared) {
                double sigma_over_r_squared =
                    sigma_squared / distance_squared;

                double sigma_over_r_6 =
                    sigma_over_r_squared
                    * sigma_over_r_squared
                    * sigma_over_r_squared;

                double sigma_over_r_12 =
                    sigma_over_r_6
                    * sigma_over_r_6;

                double pair_energy =
                    4.0
                    * epsilon
                    * (
                        sigma_over_r_12
                        - sigma_over_r_6
                    );

                total_energy += pair_energy;
            }
        }
    }

    return total_energy;
}