/* 
 * File:   Hermite.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:45 PM
 */

#ifndef HERMITE_H
#define	HERMITE_H

#include <armadillo> 
//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>

using namespace std;
using namespace arma;

/** @brief Class that contains the first Hermite polynomials, needed for the 
 * single-particle orbitals
    @author sarahrei
    @date April 2012
 */
class Hermite{
     
public:
    
    /**
     * First Hermite polynomial
     * @param x
     * @return value at x
     */
    double H_0(double x);
    
    /**
     * Second Hermite polynomial
     * @param x
     * @return value at x
     */
    double H_1(double x);
    
    /**
     * Third Hermite polynomial
     * @param x
     * @return value at x
     */
    double H_2(double x);
};



#endif	/* HERMITE_H */

