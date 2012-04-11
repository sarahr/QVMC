
#include "qvmc.h"

////////////////////////////////////////////////////////////////////////////////
//                            class Radial

/** @brief Class that computes all radial positions and distances between 
 * particles and saves them for several usages
    @author sarahrei
    @date 11 April  2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param num_part - number of particles
 * @param dim - dimension
 */
Radial::Radial(int num_part, int dim) {
    this->num_part = num_part;
    this->dim = dim;

}

/**
 * This function computes the radius squared r^2 of each particle in the 
 * vector r = (r_0^2, r_1^2,...) and the distances r_ij between all 
 * particles in the (upper-triangular!) matrix r_int
 * @param R - matrix with the Cartesian coordinates of all particles
 */
void Radial::update(mat& R) {

    int i, j, k;

    current = R;
    r.zeros(num_part);
    r_int.zeros(num_part, num_part);

    // Radial part
    for (j = 0; j < dim; j++) {
        for (i = 0; i < num_part; i++) {
            r(i) += R(i, j) * R(i, j);
        }
    }

    // Interaction matrix 

    for (k = 0; k < dim; k++) {
        for (j = 1; j < num_part; j++) {
            for (i = 0; i < j; i++) {
                r_int(i, j) += (R(j, k) - R(i, k))*(R(j, k) - R(i, k));
            }
        }
    }

    for (j = 1; j < num_part; j++) {
        for (i = 0; i < j; i++) {
            r_int(i, j) = sqrt(r_int(i, j));
        }
    }

    return;
}

/**
 *  Same computations as update(), but optimized if only particle "p" has 
 * been moved
 * @param p - particle that has been moved
 * @param R - matrix with the Cartesian coordinates of all particles
 */
void Radial::update_p(int p, mat& R) {

    int k, i;

    for (k = 0; k < dim; k++) {
        current(p, k) = R(p, k);
    }

    // Radial part
    r(p) = 0;
    for (k = 0; k < dim; k++) {
        r(p) += R(p, k) * R(p, k);
    }

    // Interaction matrix

    for (i = 0; i < p; i++) {

        r_int(i, p) = 0;

        for (k = 0; k < dim; k++) {
            r_int(i, p) += (R(p, k) - R(i, k))*(R(p, k) - R(i, k));
        }

        r_int(i, p) = sqrt(r_int(i, p));
    }


    for (i = p + 1; i < num_part; i++) {

        r_int(p, i) = 0;

        for (k = 0; k < dim; k++) {
            r_int(p, i) += (R(i, k) - R(p, k))*(R(i, k) - R(p, k));
        }

        r_int(p, i) = sqrt(r_int(p, i));
    }

    return;

}