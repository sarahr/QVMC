/* 
 * File:   QForce.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:46 PM
 */

#ifndef QFORCE_H
#define	QFORCE_H

//#include <armadillo> 
#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;
 
/** @brief Class for the quantum force
    @author sarahrei
    @date 11 April  2012
 */
class QForce{
    
private:
    Wavefunction* Psi;
    int dim, numpart;
    bool interaction;
    
public:
    
    /**
     * Constructor
     * @param Psi - Pointer to wavefunction
     */
    QForce(Wavefunction* Psi);
    
    /**
     * Compute the new quantum force after particle "p" has been moved
     * @param p - particle that has been moved
     * @param R - matrix with current position
     * @param Pos - pointer to object of class Radial with current position
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @param test_inv - inverse of Slater determinant with trial position
     */
    void new_qf(int p, mat& R, Radial* Pos, double alpha, double beta, mat& test_inv);
    
    mat qf_old;
    mat qf_new;
    
  
};


#endif	/* QFORCE_H */

