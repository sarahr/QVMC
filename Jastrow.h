/* 
 * File:   Jastrow.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:46 PM
 */

#ifndef JASTROW_H
#define	JASTROW_H

#include <armadillo> 
//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;
 
/** @brief Class for the Jastrow factor of the wavefunction 
    @author sarahrei
    @date 11 April 2012
 */
class Jastrow{
private:
    int numpart, dim;
    mat a; 
    vec grad;
   
public:
    mat g_ij; 
    vec g_new;
    Jastrow(){};
    
    /**
     * Constructor
     * @param numpart - number of particles
     * @param dim - dimension
     */
    Jastrow(int numpart,int dim); 
    
    /**
     * Compute value of Jastrow function
     * @param Pos - Pointer to object of class Radial with current position
     * @param beta - variational parameter
     * @return value of Jastrow function
     */
    double value(Radial* Pos, double beta);
    
    /**
     * Compute ratio between new and old Jastrow function after particle p has
     * been moved
     * @param Pos - pointer to object of class Radial with current position
     * @param beta - variational parameter
     * @param p - particle that has been moved
     * @return  ratio between new and old Jastrow factor
     */
    double ratio(Radial* Pos, double beta, int p);
    
    /**
     * Compute the gradient of the Jastrow function with respect to particle p
     * @param Pos - pointer to object of class Radial with current position
     * @param p- particle that has been moved
     * @param beta - variational parameter
     * @return gradient of the Jastrow function with respect to particle p
     */
    vec gradient(Radial* Pos, int p, double beta);
    
    /**
     * Computes a term needed repeatedly in the function gradient().
     * Note: "p" smaller than "i" !              
     * @param Pos - Pointer to object of class Radial with current position
     * @param p - particle that has been moved
     * @param i - index of particle, given by function gradient()
     * @param beta - variational parameter
     * @return term for function gradient()
     */
    vec grad_term(Radial* Pos, int p, int i, double beta);
       
    /**
     * Compute the Laplacian of the Jastrow function with respect to particle p.
     * @param Pos - Pointer to object of class Radial with current position
     * @param p - particle that has been moved
     * @param beta - variational parameter
     * @return Laplacian
     */
    double laplace(Radial* Pos, int p, double beta);
    
    /**
     * Computes a term needed repeatedly in the function laplace().
     * Note: "p" smaller than "i" !
     * @param Pos - Pointer to object of class Radial with current position
     * @param p - particle that has been moved
     * @param i - index of particle, given by function gradient()
     * @param beta - variational parameter
     * @return term for function laplace()
     */
    double laplace_term(Radial* Pos, int p, int i, double beta);
};


#endif	/* JASTROW_H */

