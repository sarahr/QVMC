
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

/**
 * Derivative of the first Hermite polynomial
 * @param x
 * @return derivative with respect to x
 */
double Hermite::H_0_deriv(double x) {
    return 0.0;
}

/**
 * Derivative of the second Hermite polynomial
 * @param x
 * @return derivative with respect to x
 */
double Hermite::H_1_deriv(double x) {
    return 2;
}

/**
 * Derivative of the third Hermite polynomial
 * @param x
 * @return derivative with respect to x
 */
double Hermite::H_2_deriv(double x) {
    return 8 * x;
}