/* 
 * File:   Wavefunction.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:47 PM
 */

#ifndef WAVEFUNCTION_H
#define	WAVEFUNCTION_H

#include <armadillo> 
//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;
 

/** @brief Class for the trial wavefunctions.
 * Unifies Slater part, exponential function and Jastrow function
    @author sarahrei
    @date 11 April  2012
 */
class Wavefunction{ 
private:
    int dim;
    int numpart; 
    double cur_val;
    bool correl;
           
public:
    
    /**
     * Constructor
     * @param numpart - number of particles
     * @param dim - dimension
     * @param omega - oscillator frequency
     * @param jastrow - 0: without Jastrow factor, 1: with Jastrow factor
     */
    Wavefunction(int numpart, int dim, double omega, int jastrow);
    
    /**
     * Compute value of the wavefunction
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return value of the wavefunction
     */
    double value(mat& R, double alpha, double beta);
    
    /**
     * Compute the ratio between new and old wavefunction after particle p
     * has been moved. The main computations are done in the subclasses, here 
     * everything is set together.
     * @param R_new - matrix containing the new Cartesian coordinates of all particles
     * @param p - particle that has been moved
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return ratio between new and old wavefunction after particle p
     * has been moved
     */
    double ratio(mat& R_new, int p, double alpha, double beta);
    
    /**
     * @return current value of the wavefunction
     */
    double get_curval() const;
    
    /**
     * Returns whether or not the correlation term is considered
     * @return whether or not the correlation term is considered
     */
    bool get_int()const;
    
    /**
     * @return dimension
     */
    int get_dim() const;
    
    /**
     * @return number of particles
     */
    int get_numpart() const;
    
    Slater* SlaterPsi;
    ExpFactor* ExpFactorPsi;
    Jastrow* JastrowPsi;
    Radial* Pos, *Pos_tr;
   
};


#endif	/* WAVEFUNCTION_H */

