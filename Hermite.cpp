
#include "qvmc.h"

////////////////////////////////////////////////////////////////////////////////
//                          class Hermite,
//                         
/** @brief Class that contains the first Hermite polynomials, needed for the 
 * single-particle orbitals
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * First Hermite polynomial
 * @param x
 * @return value at x
 */
double Hermite::H_0(double x) {
    return 1.0;
}

/**
 * Second Hermite polynomial
 * @param x
 * @return value at x
 */
double Hermite::H_1(double x) {
    return 2 * x;
}

/**
 * Third Hermite polynomial
 * @param x
 * @return value at x
 */
double Hermite::H_2(double x) {
    return 4 * x * x - 2;
}
