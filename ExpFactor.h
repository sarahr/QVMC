/* 
 * File:   ExpFactor.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:45 PM
 */

#ifndef EXPFACTOR_H
#define	EXPFACTOR_H

#include <armadillo> 
//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;

/** @brief Class that handles the exponential factor of the wavefunction
    @author sarahrei
    @date 11 April 2012
 */
class ExpFactor {
private:
    int numpart;
    double omega;
    double t_value;

public:

    ExpFactor();

    /**
     * Constructor
     * @param numpart - number of particles
     * @param omega - osc. frequency omega
     */
    ExpFactor(int numpart, double omega);

    /**Compute the value of the exponential factor of the wavefunction
        @param Pos - Radial-object of the current position
        @param alpha - variational parameter
        @return value of the exponential factor of the wavefunction
     */
    double value(Radial* Pos, double alpha);

    /**Compute the ratio between new and old exponential factor after particle p
   has been moved.
        @param Pos - Radial-object of the previous position
        @param Pos_tr - Radial-object of the new position
        @param alpha - variational parameter
        @param p - moved particle
        @return ratio between new and old exponential factor
     */
    double ratio(Radial* Pos, Radial* Pos_tr, double alpha, int p);
};


#endif	/* EXPFACTOR_H */
