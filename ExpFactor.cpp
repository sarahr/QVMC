
#include "qvmc.h" 

////////////////////////////////////////////////////////////////////////////////
//                             class ExpFactor

/** @brief Class that handles the exponential factor of the wavefunction
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 * @param omega - osc. frequency omega
 */
ExpFactor::ExpFactor(int numpart, double omega) {
    this->numpart = numpart;
    this->omega = omega;

}

/**Compute the value of the exponential factor of the wavefunction
        @param Pos - Radial-object of the current position
        @param alpha - variational parameter
        @return value
 */
double ExpFactor::value(Radial* Pos, double alpha) {

    double term = 0.0;

    for (int i = 0; i < numpart; i++) {
        term += Pos->r(i);
    }

    term *= -0.5 * alpha*omega;
    t_value = exp(term);

    return t_value;
}

/**Compute the ratio between new and old exponential factor after particle p
   has been moved.
        @param Pos - Radial-object of the previous position
        @param Pos_tr - Radial-object of the new position
        @param alpha - variational parameter
        @param p - moved particle
        @return ratio
 */
double ExpFactor::ratio(Radial* Pos, Radial* Pos_tr, double alpha, int p) {

    double term;

    term = Pos_tr->r(p) - Pos->r(p);
    term *= -0.5 * alpha*omega;
    term = exp(term);

    return term;
}
