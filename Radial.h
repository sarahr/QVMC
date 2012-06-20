/* 
 * File:   Radial.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:47 PM
 */

#ifndef RADIAL_H
#define	RADIAL_H

//#include <armadillo> 
#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;
 
/** @brief Class that computes all radial positions and distances between 
 * particles and saves them for several usages
    @author sarahrei
    @date 11 April  2012
 */
class Radial { 
public:

    Radial() {};
    
    /**
     * Constructor
     * @param num_part - number of particles
     * @param dim - dimension
     */
    Radial(int num_part, int dim);
    vec r;
    mat r_int, current;
    int dim, num_part;

    /**
     * This function computes the radius squared r^2 of each particle in the 
     * vector r = (r_0^2, r_1^2,...) and the distances r_ij between all 
     * particles in the (upper-triangular!) matrix r_int
     * @param R - matrix with the Cartesian coordinates of all particles
     */
    void update(mat& R);
    
    /**
     *  Same computations as update(), but optimized if only particle "p" has 
     * been moved
     * @param p - particle that has been moved
     * @param R - matrix with the Cartesian coordinates of all particles
     */
    void update_p(int p, mat& R);
}; 


#endif	/* RADIAL_H */

