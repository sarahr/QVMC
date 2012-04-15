
#include "qvmc.h" 


////////////////////////////////////////////////////////////////////////////////
//                          class Wavefunction 

/** @brief Class for the trial wavefunctions
    @author sarahrei
    @date 11 April  2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 * @param dim - dimension
 * @param omega - oscillator frequency
 * @param jastrow - 0: without Jastrow factor, 1: with Jastrow factor
 */
Wavefunction::Wavefunction(int numpart, int dim, double omega, int jastrow) {

    this->dim = dim;
    this->numpart = numpart;

    if (jastrow == 1) correl = true; // with or without correlation term
    else correl = false;

    SlaterPsi = new Slater(numpart, omega, dim);
    ExpFactorPsi = new ExpFactor(numpart, omega);
    JastrowPsi = new Jastrow(numpart, dim);
    Pos = new Radial(numpart, dim);
    Pos_tr = new Radial(numpart, dim);

}

/**
 * Compute value of the wavefunction
 * @param R - matrix containing the Caresian coordinates of all particles
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return value
 */
double Wavefunction::value(mat& R, double alpha, double beta) {

    cur_val = 1.0;
    Pos_tr->update(R);

    // Get the different contributions
    cur_val *= SlaterPsi->value(R, alpha);
    cur_val *= ExpFactorPsi->value(Pos_tr, alpha);

    // Jastrow factor only if correlation
    if (correl) cur_val *= JastrowPsi->value(Pos_tr, beta);

    return cur_val;

}

/**
 * Compute the ratio between new and old wavefunction after particle p
 * has been moved. The main computations are done in the subclasses, here 
 * everything is set together.
 * @param R_new - matrix containing the new Cartesian coordinates of all particles
 * @param p - particle that has been moved
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return ratio
 */
double Wavefunction::ratio(mat& R_tr, int p, double alpha, double beta) {

    double r = 1.0;

    r *= SlaterPsi->ratio(R_tr, p, alpha);
    r *= ExpFactorPsi->ratio(Pos, Pos_tr, alpha, p);
    if (correl) r *= JastrowPsi->ratio(Pos_tr, beta, p);

    SlaterPsi->update_inverse(p);

    return r;

}

/**
 * Returns whether or not the correlation term is considered
 * @return correlation
 */
bool Wavefunction::get_int() const {
    return correl;
}

/**
 * @return current value of the wavefunction
 */
double Wavefunction::get_curval() const {
    return cur_val;
}

/**
 * @return dimension
 */
int Wavefunction::get_dim() const {
    return dim;

}

/**
 * @return number of particles
 */
int Wavefunction::get_numpart() const {
    return numpart;
}


