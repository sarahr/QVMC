
#include "qvmc.h" 

////////////////////////////////////////////////////////////////////////////////
//                           class Slater

/** @brief Class for the slater determinant part of the wavefunction (single-
 * particle orbitals without exponential factor
    @author sarahrei
    @date 11 April  2012
 */
//////////////////////////////////////////////////////////////////////////////// 

/**
 * Constructor
 * @param numpart - number of particles 
 * @param omega - osc. frequency omega
 * @param dim - dimension
 */
Slater::Slater(int numpart, double omega, int dim) {
    this->numpart = numpart;
    this->omega = omega;
    this->dim = dim;

    sqom = sqrt(omega);
    n2 = numpart / 2;

    // Initialization, setting size of vectors/matrices
    slater = zeros(n2, numpart);
    slat_inv = ones(n2, numpart);
    inv_back = zeros(n2, numpart);
    detv = zeros(2); // spin up/spin down determinant
    grad = zeros(dim);
    sinp_new = zeros(n2);

    // Set up the ordering of the single-particle states
    part_code.set_size(6, 2);
    part_code << 0 << 0 << endr << 1 << 0 << endr << 0 << 1 << endr << 2 << 0 <<
            endr << 1 << 1 << endr << 0 << 2 << endr;

}


/*
 * Array of functions containing the Hermite polynomials, needed for the 
 * single-particle functions
 */
const Slater::fptr Slater::herm_table[] = {
    &Hermite::H_0, &Hermite::H_1, &Hermite::H_2
};

/**
 * Compute the value of the Slater determinant
 * @param R - matrix containing the Cartesian coordinates of all particles
 * @param alpha - variational parameter
 * @return value
 */
double Slater::value(mat& R, double alpha) {

    double sqal = sqrt(alpha);

    for (int p = 0; p < n2; p++) { // Loop over all particles
        for (int i = 0; i < n2; i++) { // Loop over all single-particle states

            // Spin up
            slater(i, p) = (BasSin->*herm_table[part_code(i, 0)])
                    (sqal * sqom * R(p, 0))*(BasSin->*herm_table[part_code(i, 1)])
                    (sqal * sqom * R(p, 1));

            // Spin down
            slater(i, p + n2) = (BasSin->*herm_table[part_code(i, 0)])
                    (sqal * sqom * R(p + n2, 0))*(BasSin->*herm_table[part_code(i, 1)])
                    (sqal * sqom * R(p + n2, 1));
        }
    }

    detv(0) = det(slater.submat(0, 0, n2 - 1, n2 - 1)); // Spin up
    detv(1) = det(slater.submat(0, n2, n2 - 1, numpart - 1)); // Spin down

    return detv(0) * detv(1);

}

/**
 * Compute the gradient of the Slater determinant with respect to particle "p"
 * Note: Gradient is of complete SD, including exponential factor!
 * @param R - matrix containing the Cartesian coordinates of all particles
 * @param p - particle that has been moved
 * @param alpha - variational parameter
 * @param inv - Slater inverse
 * @return gradient
 */
vec Slater::gradient(mat& R, int p, double alpha, mat& inv) {

    grad.zeros();

    for (int i = 0; i < n2; i++) {
        grad += orb_grad(R, p, i, alpha) * inv(i, p);
    }

    return grad;
}

/**
 * Computes the gradient of the single-particle orbitals, needed by the 
 * function gradient()
 * @param R - matrix containing the Cartesian coordinates of all particles
 * @param p - particle that has been moved, supplies coordinates
 * @param qq  - supplies quantum numbers
 * @param alpha - variational parameter
 * @return gradient of the single-particle orbitals
 */
vec Slater::orb_grad(mat& R, int p, int qq, double alpha) {

    vec gra_ob(dim);
    double sqal = sqrt(alpha);
    double f = 2 * sqom*sqal;

    for (int i = 0; i < dim; i++) {
        gra_ob(i) = f * part_code(qq, i)*(BasSin->*herm_table[abs(part_code(qq, i) - 1)])
                (sqal * sqom * R(p, i)); // abs to prevent negative indices, in this
        // case the product gets zero anyway

        gra_ob(i) /= (BasSin->*herm_table[part_code(qq, i)])(sqal * sqom * R(p, i));
        gra_ob(i) -= R(p, i) * omega*alpha;
    }

    gra_ob *= (BasSin->*herm_table[part_code(qq, 0)])(sqal * sqom * R(p, 0))*
            (BasSin->*herm_table[part_code(qq, 1)])(sqal * sqom * R(p, 1));

    return gra_ob;

}

/**
 * Compute the Laplacian of the Slater determinant with respect to particle "p"
 * Note: Laplacian is of complete SD, including exponential factor!
 * @param Pos - pointer to object of class Radial with current position
 * @param p - particle that has been moved, supplies coordinates
 * @param alpha - variational parameter
 * @return Laplacian
 */
double Slater::laplace(Radial* Pos, int p, double alpha) {

    double lapl = 0.0;

    for (int i = 0; i < n2; i++) {
        lapl += orb_laplace(Pos, p, i, alpha) * slat_inv(i, p);
    }

    return lapl;

}

/**
 * Computes the Laplacian of the single-particle orbitals, needed by the 
 * function laplace()
 * @param Pos - pointer to object of class Radial with current position
 * @param p - particle that has been moved, supplies coordinates
 * @param q - determines single-particle orbital
 * @param alpha - variational parameter
 * @return Laplacian of the single-particle orbitals
 */
double Slater::orb_laplace(Radial* Pos, int p, int qq, double alpha) {

    double l;
    double sqal = sqrt(alpha);
    int nx = part_code(qq, 0);
    int ny = part_code(qq, 1);

    l = 4 * nx * (nx - 1)*(BasSin->*herm_table[abs(part_code(qq, 0) - 2)])
            (sqal * sqom * Pos->current(p, 0)) / ((BasSin->*herm_table[part_code(qq, 0)])
            (sqal * sqom * Pos->current(p, 0)));
    l += 4 * ny * (ny - 1)*(BasSin->*herm_table[abs(part_code(qq, 1) - 2)])
            (sqal * sqom * Pos->current(p, 1)) / ((BasSin->*herm_table[part_code(qq, 1)])
            (sqal * sqom * Pos->current(p, 1)));

    l += omega * alpha * Pos->r(p);

    l -= 4 * nx * Pos->current(p, 0) * sqal * sqom * (BasSin->*herm_table[abs(part_code(qq, 0) - 1)])
            (sqal * sqom * Pos->current(p, 0)) / ((BasSin->*herm_table[part_code(qq, 0)])
            (sqal * sqom * Pos->current(p, 0)));
    l -= 4 * ny * Pos->current(p, 1) * sqal * sqom * (BasSin->*herm_table[abs(part_code(qq, 1) - 1)])
            (sqal * sqom * Pos->current(p, 1)) / ((BasSin->*herm_table[part_code(qq, 1)])
            (sqal * sqom * Pos->current(p, 1)));
    l -= 2;
    l *= alpha*omega;

    l *= (BasSin->*herm_table[part_code(qq, 0)])(sqal * sqom * Pos->current(p, 0))
            *(BasSin->*herm_table[part_code(qq, 1)])(sqal * sqom * Pos->current(p, 1));

    return l;

}

/**
 * Compute the complete inverse of the spin-up or spin-down part of the 
 * Slater determinant
 * @param p - determines which spin part is calculated
 * @param R - matrix containing the Cartesian coordinates of all particles
 * @param alpha - variational parameter
 * @return Slater inverse 
 */
mat Slater::inverse(int p, mat& R, double alpha) {

    int access = p / n2; // Integer division to access correct spin
    mat test_inv = ones(n2, numpart);
    value(R, alpha);

    test_inv.submat(0, access*n2, n2 - 1, n2 - 1 + access * n2) = inv(slater.submat(0, access*n2, n2 - 1, n2 - 1 + access * n2));

    return test_inv;

}

/**
 * Update the inverse of the Slater determinant
 * Optimized version for movement of one particle at a time
 * @param p - particle that has been moved
 */
void Slater::update_inverse(int p) {

    double s;
    int access = p / n2; // to access right spin


    // First backup in case move is rejected
    for (int j = access * n2; j < (access + 1) * n2; j++) {
        for (int k = 0; k < n2; k++) {
            inv_back(k, j) = slat_inv(k, j);
        }
    }


    // First update j \neq p
    for (int j = access * n2; j < p; j++) {
        s = 0.0;

        for (int k = 0; k < n2; k++) {
            s += sinp_new(k) * slat_inv(k, j);
        }

        for (int k = 0; k < n2; k++) {
            slat_inv(k, j) -= (s / cur_rat) * slat_inv(k, p);
        }
    }

    for (int j = p + 1; j < (access + 1) * n2; j++) {
        s = 0.0;

        for (int k = 0; k < n2; k++) {
            s += sinp_new(k) * slat_inv(k, j);
        }

        for (int k = 0; k < n2; k++) {
            slat_inv(k, j) -= (s / cur_rat) * slat_inv(k, p);
        }
    }

    // now for j= p
    for (int k = 0; k < n2; k++) {
        slat_inv(k, p) /= cur_rat;
    }

    return;

}

/**
 * Compute the ratio between new and old Slater determinant
 * @param R_new - matrix containing Cartesian coordiantes of all particles
 * @param p  - determines which of the two determinants (spin up or down)
 * is considered
 * @param alpha - variational parameter
 * @return ratio 
 */
double Slater::ratio(mat& R_new, int p, double alpha) {

    cur_rat = 0.0;
    double sqal = sqrt(alpha);

    for (int j = 0; j < n2; j++) {
        sinp_new(j) = (BasSin->*herm_table[part_code(j, 0)])(sqal * sqom * R_new(p, 0))
                *(BasSin->*herm_table[part_code(j, 1)])(sqal * sqom * R_new(p, 1));
    }

    for (int j = 0; j < n2; j++) {
        cur_rat += sinp_new(j) * slat_inv(j, p);
    }

    return cur_rat;

}

/**
 * Resets the Slater inverse in case the move has been rejected
 * @param p - particle that has been moved
 */
void Slater::reset(int p) {

    int access = p / n2; // access correct spin

    for (int j = access * n2; j < (access + 1) * n2; j++) {
        for (int k = 0; k < n2; k++) {
            slat_inv(k, j) = inv_back(k, j);
        }
    }

    return;

}
