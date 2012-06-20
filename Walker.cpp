
#include "Walker.h" 
#include "lib.h"
#include <time.h>


////////////////////////////////////////////////////////////////////////////////
//                             class Walker

/** @brief Class for handling the walkers in the DMC algorithm
    @author sarahrei
    @date 15 May 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 * @param dim - dimension
 * @param omega - oscillator frequency
 * @param jastrow - 0: without, 1: with Jastrow factor
 * @param H - Hamiltonian of the system
 * @param alpha - first variational parameter of trial wavefunction
 * @param beta - second variational parameter of trial wavefunction
 */
Walker::Walker(int numpart, int dim, double omega, int jastrow,
        Hamiltonian* H, double alpha, double beta) {

    this->numpart = numpart;
    this->dim = dim;
    this->alpha = alpha;
    this->beta = beta;

    Trial = new Wavefunction(numpart, dim, omega, jastrow);
    QFo = new QForce(Trial);

    this->H = H;
    interaction = H->get_interaction();

    R_cur.set_size(numpart, dim);
    R_tr.set_size(numpart, dim);

}

/**
 * Destructor
 */
Walker::~Walker() {

    delete Trial;
    delete QFo;
}

/**
 * Initialization of walker
 * @param del_t - time step
 */
void Walker::initWalker(double del_t) {

    this->del_t = del_t;

    // Initialize the position
    R_cur = randu(numpart, dim);
    R_tr = R_cur;
    Trial->Pos->update(R_cur);
    Trial->Pos_tr->update(R_tr);

    // Initialize Slater determinant and Jastrow factor
    Trial->SlaterPsi->value(R_cur, alpha);
    Trial->JastrowPsi->value(Trial->Pos, beta);

    // Initialize Slater inverse
    int n2 = numpart / 2;
    for (int i = 0; i < 2; i++) { // For spin up and down
        Trial->SlaterPsi->slat_inv.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2)
                = inv(Trial->SlaterPsi->slater.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2));
    }
    trial_inv = Trial->SlaterPsi->slat_inv;

    for (int p = 0; p < numpart; p++) {
        QFo->new_qf(p, R_cur, Trial->Pos, alpha, beta, trial_inv);
    }

    QFo->qf_old = QFo->qf_new;

    return;

}

/**
 * Compute the local energy
 * @return local energy
 */
double Walker::E_local() {

    double energy;

    energy = H->H_0(Trial, alpha, beta);
    if (interaction) energy += H->H_1(Trial); // interaction part optional

    return energy;

}

/**
 * If the move is not accepted, reset position and Slater inverse
 * @param p - particle that has been moved
 */
void Walker::not_accept(int p) {

    for (int k = 0; k < dim; k++) {
        R_tr(p, k) = R_cur(p, k);
    }

    Trial->Pos_tr->update_p(p, R_cur);
    Trial->SlaterPsi->reset(p);

    return;

}

/**
 * Compute trial position according to the Metropolis-Hastings algorithm
 * @param p - particle that has been moved 
 
 */
void Walker::trial_pos(int p) {

    double chi;

    for (int j = 0; j < dim; j++) {
        chi = DRanNormalZigVec();
        R_tr(p, j) = R_cur(p, j) + 0.5 * QFo->qf_old(p, j) * del_t + chi * sqrt(del_t);

    }

    Trial->Pos_tr->update_p(p, R_tr);
    return;

}

/**
 * Compute acceptance ratio according to the Metropolis-Hastings algorithm
 * @param p - particle that has been moved 
 * @return acceptance ratio 
 */
double Walker::ratio(int p) {

    double r = 1.0, green = 0.0;
    vec ef1(dim), ef2(dim);

    // Ratio between SDs to check if node has been crossed
    wf_R = Trial->SlaterPsi->ratio(R_tr, p, alpha);
    r *= wf_R;

    if (Trial->correl) r *= Trial->JastrowPsi->ratio(Trial->Pos_tr, beta, p);
    r *= r;

    Trial->SlaterPsi->update_inverse(p);

    QFo->new_qf(p, R_tr, Trial->Pos_tr, alpha, beta, Trial->SlaterPsi->slat_inv);

    for (int i = 0; i < numpart; i++) {
        for (int k = 0; k < dim; k++) {

            if (i != p) {
                green += 0.25 * del_t * 0.5 *
                        (QFo->qf_old(i, k) + QFo->qf_new(i, k)) *
                        (QFo->qf_old(i, k) - QFo->qf_new(i, k));
            } else {

                green += 0.5 * (QFo->qf_old(i, k) + QFo->qf_new(i, k)) *
                        (0.5 * 0.5 * del_t * (QFo->qf_old(i, k) - QFo->qf_new(i, k)) +
                        R_cur(i, k) - R_tr(i, k));
            }
        }
    }

    green = exp(green);
    r *= green;

    return r;

}

/**
 * Performs all updates if a move has been accepted
 * @param p - particle that has been moved 
 */
void Walker::accept(int p) {

    for (int k = 0; k < dim; k++) {
        R_cur(p, k) = R_tr(p, k); // Update position
    }

    QFo->qf_old = QFo->qf_new;

    Trial->Pos->update_p(p, R_cur);

    // Update Jastrow factor
    if (interaction) {
        for (int i = 0; i < p; i++) {
            Trial->JastrowPsi->g_ij(i, p) = Trial->JastrowPsi->g_new(i);
        }

        for (int i = p + 1; i < numpart; i++) {
            Trial->JastrowPsi->g_ij(p, i) = Trial->JastrowPsi->g_new(i);
        }
    }

    return;

}

/**
 * Check if a node has been crossed
 * @return "true" if node crossed, otherwise "false"
 */
bool Walker::nodeCrossed() {

    if (wf_R <= 0)
        return true;
    else
        return false;

}


