
#include "qvmc.h"  


//////////////////////////////////////////////////////////////////////////////// 
//                           class Jastrow 

/** @brief Class for the Jastrow factor of the wavefunction 
    @author sarahrei
    @date 11 April  2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 * @param dim - dimension
 */
Jastrow::Jastrow(int numpart, int dim) {
    this->numpart = numpart;
    this->dim = dim;

    int i, j;
    int n2 = numpart / 2;

    a.zeros(numpart, numpart);
    g_ij.zeros(numpart, numpart);
    g_new.zeros(numpart);
    grad.zeros(dim);

    for (j = 0; j < numpart; j++) { // Diagonal elements
        a(j, j) = 1. / 3;
    }

    for (j = 0; j < n2; j++) { // Equal spin
        for (i = 0; i < j; i++) {
            a(i, j) = 1. / 3;
            a(i + n2, j + n2) = 1. / 3;
        }
    }

    for (j = n2; j < numpart; j++) { // Opposite spin
        for (i = 0; i < n2; i++) {
            a(i, j) = 1.0;
        }
    }

}

/**
 * Compute value of Jastrow function
 * @param Pos - Pointer to object of class Radial with current position
 * @param beta - variational parameter
 * @return value
 */
double Jastrow::value(Radial* Pos, double beta) {

    double ret;
    double sum = 0.0;

    for (int j = 1; j < numpart; j++) {
        for (int i = 0; i < j; i++) {

            g_ij(i, j) = a(i, j) * Pos->r_int(i, j) / (1 + beta * Pos->r_int(i, j));
            sum += g_ij(i, j);
        }
    }

    ret = exp(sum);

    return ret;
}

/**
 * Compute ratio between new and old Jastrow function after particle p has
 * been moved
 * @param Pos - Pointer to object of class Radial with current position
 * @param beta - variational parameter
 * @param p - particle that has been moved
 * @return  ratio between new and old Jastrow factor
 */
double Jastrow::ratio(Radial* Pos, double beta, int p) {

    double rj = 0.0;

    for (int i = 0; i < p; i++) {

        g_new(i) = a(i, p) * Pos->r_int(i, p) / (1 + beta * Pos->r_int(i, p));
        rj += g_new(i) - g_ij(i, p);
    }

    for (int i = p + 1; i < numpart; i++) {
        g_new(i) = a(p, i) * Pos->r_int(p, i) / (1 + beta * Pos->r_int(p, i));
        rj += g_new(i) - g_ij(p, i);

    }

    rj = exp(rj);

    return rj;

}

/**
 * Compute the gradient of the Jastrow function with respect to particle p
 * @param Pos - Pointer to object of class Radial with current position
 * @param p- particle that has been moved
 * @param beta - variational parameter
 * @return gradient
 */
vec Jastrow::gradient(Radial* Pos, int p, double beta) {

    vec grad = zeros(dim);

    for (int k = 0; k < p; k++) {
        grad += grad_term(Pos, p, k, beta);
    }

    for (int k = p + 1; k < numpart; k++) {
        grad += grad_term(Pos, p, k, beta);
    }

    return grad;
}

/**
 * Computes a term needed repeatedly in the function gradient().
 * Note: p<i!
 * @param Pos - Pointer to object of class Radial with current position
 * @param p - particle that has been moved
 * @param i - index of particle, given by function gradient()
 * @param beta - variational parameter
 * @return term for function gradient()
 */
vec Jastrow::grad_term(Radial* Pos, int p, int i, double beta) {

    vec ret(dim);

    for (int j = 0; j < dim; j++) {
        ret(j) = Pos->current(p, j) - Pos->current(i, j);
    }

    if (p > i) {
        int tmp = p;
        p = i;
        i = tmp;
    }

    double deno = 1 + beta * Pos->r_int(p, i);
    ret *= a(p, i) / (Pos->r_int(p, i) * deno * deno);

    return ret;

}

/**
 * Compute the Laplacian of the Jastrow function with respect to particle p.
 * @param Pos - Pointer to object of class Radial with current position
 * @param p - particle that has been moved
 * @param beta - variational parameter
 * @return Laplacian
 */
double Jastrow::laplace(Radial* Pos, int p, double beta) {

    double ret = 0.0;

    ret += dot(gradient(Pos, p, beta), gradient(Pos, p, beta));

    for (int k = 0; k < p; k++) {
        ret += laplace_term(Pos, k, p, beta);
    }

    for (int k = p + 1; k < numpart; k++) {
        ret += laplace_term(Pos, p, k, beta);
    }

    return ret;
}

/**
 * Computes a term needed repeatedly in the function laplace().
 * Note: p<i
 * @param Pos - Pointer to object of class Radial with current position
 * @param p - particle that has been moved
 * @param i - index of particle, given by function gradient()
 * @param beta - variational parameter
 * @return term for function laplace()
 */
double Jastrow::laplace_term(Radial* Pos, int p, int i, double beta) {

    double ret = 1.0;
    double r_ip = Pos->r_int(p, i);

    double deno = 1 + beta*r_ip;
    ret *= (1 - beta * r_ip);
    ret *= a(p, i) / (r_ip * deno * deno * deno);

    return ret;

}
