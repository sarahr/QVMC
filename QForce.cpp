#include "qvmc.h"  

////////////////////////////////////////////////////////////////////////////////
//                           class QForce

/** @brief Class for the quantum force
    @author sarahrei
    @date 15 May 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param Psi - Pointer to wavefunction
 */
QForce::QForce(Wavefunction* Psi) {
    this->Psi = Psi;
    this->dim = Psi->get_dim();
    this->numpart = Psi->get_numpart();

    qf_old = zeros(numpart, dim);
    qf_new = zeros(numpart, dim);

    interaction = Psi->get_int();

}

/**
 * Compute the new quantum force after particle "p" has been moved
 * @param p - particle that has been moved
 * @param R - matrix with current position
 * @param Pos - pointer to object of class Radial with current position
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @param trial_inv - inverse of Slater determinant with trial position
 */
void QForce::new_qf(int p, mat& R, Radial* Pos, double alpha, double beta,
        mat& trial_inv) {

    int n2 = numpart / 2;

    for (int pp = 0; pp < numpart; pp++) {
        for (int i = 0; i < dim; i++) {

            qf_new(pp, i) = (Psi->SlaterPsi->gradient(R, pp, alpha, trial_inv))(i);
            if (interaction) qf_new(pp, i) += (Psi->JastrowPsi->gradient(Pos, pp, beta))(i);
            qf_new(pp, i) *= 2;
        }
    }

    return;

}

