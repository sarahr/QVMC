/* 
 * File:   Slater.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:47 PM
 */

#ifndef SLATER_H
#define	SLATER_H

#include <armadillo> 
//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;
 

/** @brief Class for the slater determinant part of the wavefunction (single-
 * particle orbitals without exponential factor
    @author sarahrei
    @date 11 April  2012
 */
class Slater{
private: 
    int numpart; 
    int dim;
    double omega, sqom, cur_rat; 
    mat  sd_dalpha;
    vec detv, grad, sinp_new;
    int n2;

    // Definitions for the orbital array
    typedef double (Hermite::*fptr)(double); 
    static const fptr herm_table[3];
    static const fptr deriv_table[3];

    Hermite* BasSin;
    imat part_code; // code for the single-particle basis (up to N=12 particles)
  
     
public:
    mat slater, slat_inv, inv_back;
    Slater(){};
    
    /**
     * Constructor
     * @param numpart - number of particles 
     * @param omega - osc. frequency omega
     * @param dim - dimension
     */
    Slater(int numpart, double omega, int dim);
    
    /**
     * Compute the value of the Slater determinant
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param alpha - variational parameter
     * @return value of the Slater determinant
     */
    double value(mat& R, double alpha);
    
    /**
     * Compute the ratio between new and old Slater determinant
     * @param R_new - matrix containing Cartesian coordiantes of all particles
     * @param p  - determines which of the two determinants (spin up or down)
     * is considered
     * @param alpha - variational parameter
     * @return ratio between new and old Slater determinant
     */
    double ratio(mat& R_new, int p, double alpha);
    
    /**
     * Compute the gradient of the Slater determinant with respect to particle "p"
     * Note: Gradient is of complete SD, including exponential factor!
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param p - particle that has been moved
     * @param alpha - variational parameter
     * @param inv - Slater inverse
     * @return gradient of the Slater determinant with respect to particle "p"
     */
    vec gradient(mat& R, int p, double alpha, mat& inv);
    
    /**
     * Computes the gradient of the single-particle orbitals, needed by the 
     * function gradient()
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param p - particle that has been moved, supplies coordinates
     * @param qq  - determines single-particle orbital
     * @param alpha - variational parameter
     * @return gradient of the single-particle orbitals
     */
    vec orb_grad(mat& R, int p, int qq, double alpha);
    
    /**
     * Compute the Laplacian of the Slater determinant with respect to particle "p"
     * Note: Laplacian is of complete SD, including exponential factor!
     * @param Pos - instance of "Radial" with current position
     * @param p - particle that has been moved, supplies coordinates
     * @param alpha - variational parameter
     * @return Laplacian of the Slater determinant with respect to particle "p"
     */
    double laplace(Radial* Pos, int p, double alpha);
    
    /**
     * Computes the Laplacian of the single-particle orbitals, needed by the 
     * function laplace()
     * @param Pos - instance of "Radial" with current position
     * @param p - particle that has been moved, supplies coordinates
     * @param q - determines single-particle orbital
     * @param alpha - variational parameter
     * @return Laplacian of the single-particle orbitals
     */
    double orb_laplace(Radial* Pos, int p, int q, double alpha);
    
    /**
     * Update the inverse of the Slater determinant
     * Optimized version for movement of one particle at a time
     * @param p - particle that has been moved
     */
    void update_inverse(int p);
    
    /**
     * Compute the complete inverse of the spin-up or spin-down part of the 
     * Slater determinant
     * @param p - determines which spin part is calculated
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param alpha - variational parameter
     * @return Slater inverse 
     */
    mat inverse(int p, mat& R, double alpha);
    
    /**
     * Computes the derivative of the single-particle orbitals with respect to
     * alpha, needed by the function deriv_alpha()
     * Note: derivative includes exponential factor!
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param Pos - instance of "Radial" with current position
     * @param alpha - variational parameter
     * @param q - determines single-particle orbital
     * @param p - particle that supplies coordinates
     * @return derivative with respect to alpha
     */
    double orb_deriv_alpha(mat& R, Radial* Pos, double alpha, int q, int p);
    
    /**
     * Computes the derivative of the entire Slater determinant w.r.t. alpha
     * @param R - matrix containing the Cartesian coordinates of all particles
     * @param Pos - instance of "Radial" with current position
     * @param alpha - variational parameter
     * @return derivative with respect to alpha
     */
    double deriv_alpha(mat& R, Radial* Pos, double alpha);
    
    /**
     * Resets the Slater inverse in case the move has been rejected
     * @param p - particle that has been moved
     */
    void reset(int p);
};


#endif	/* SLATER_H */

