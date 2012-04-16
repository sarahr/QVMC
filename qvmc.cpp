/*
 * Performance of the QVMC algorithm. All the contributions to the local energy
 * are written to file for further analysis by the program blocking_analyze.
 * 
 * Input arguments: 
 * argv[1]: # processors used when executing blocking.cpp
 * argv[2]: # MC cycles
 */


#include "qvmc.h"
#include "lib.h"  
#include "zignor.h" 
#include "zigrandom.h"
#include "ziggurat.hpp"
#include "normal.hpp" 

#define E_POT_KIN false // set "true" for analyzing E_pot & E_kin separately
#define MINIMIZE true // set "true" for CGM
#define DISTANCE false // set "true" if mean distance between two particles
// is to be computed



////////////////////////////////////////////////////////////////////////////////
//                             class VMC

/** @brief Class for the VMC algorithm
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param N - number of MC cycles
 * @param N_therm - number of thermalization steps
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void VMC::run_algo(int N, int N_therm, double alpha, double beta) {

    int i;
    double rat, eps;
    double del_E = 0;
    double del_Epot = 0;
    double del_Ekin = 0;


    //*************************  Thermalization  ******************************

    accepted = 0;

    initialize(alpha, beta); // Initialize the system

    for (i = 0; i < N_therm; i++) {

        for (int p = 0; p < numpart; p++) {// Loop over all particles

            // Calculate trial position
            trial_pos(p, alpha, beta);

            // Compute acceptance ratio
            rat = ratio(p, alpha, beta);

            // Check if move is accepted
            if (rat >= 1.0) { // accept if probability is greater
                accept(p, alpha, beta);
                accepted++;
            } else { // otherwise check against random number
                eps = ran3(&idum);
                if (eps < rat) {
                    accept(p, alpha, beta);
                    accepted++;
                } else
                    not_accept(p); // Do not accept
            }
        }
    }


    //********************** After thermalization *****************************


    if (N > 0) accepted = 0; // For function delta_opt()

    for (i = 0; i < N; i++) {

        for (int p = 0; p < numpart; p++) { // Loop over particles

            // Calculate trial position
            trial_pos(p, alpha, beta);

            // Compute acceptance ratio
            rat = ratio(p, alpha, beta);

            // Check if move is accepted
            if (rat >= 1.0) { // accept if probability is greater
                accept(p, alpha, beta);
                del_E = E_local(alpha, beta);
#if E_POT_KIN
                E_Pot_Kin(alpha, beta, del_Epot, del_Ekin);
#endif

#if MINIMIZE
                //part_psi(alpha, beta);
                part_psi_beta_analytic(alpha, beta);
#endif  
                accepted++;
            } else { // otherwise check against random number

                eps = ran3(&idum);
                if (eps < rat) {
                    accept(p, alpha, beta);
                    del_E = E_local(alpha, beta);
#if E_POT_KIN
                    E_Pot_Kin(alpha, beta, del_Epot, del_Ekin);
#endif

#if MINIMIZE
                    //part_psi(alpha, beta);
                    part_psi_beta_analytic(alpha, beta);
#endif   
                    accepted++;
                } else
                    not_accept(p); // Do not accept

            }

            // Updating statistics 
#if E_POT_KIN
            update_statistics(del_E, del_Epot, del_Ekin);
#else
            update_statistics(del_E);
#endif

        }
    }

    //cout << accepted/(N*numpart) << endl;

    return;

}

/**
 * Compute the local energy
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return local energy
 */
double VMC::E_local(double alpha, double beta) {

    double energy;

    energy = H->H_0(Trial, alpha, beta);
    if (interaction) energy += H->H_1(Trial); // interaction part optional

    return energy;
}

/**
 * Compute kinetic and potential energy separately
 * @param alpha - first variational parameter
 * @param beta  - second variational parameter
 * @param E_pot - variable to save the potential energy
 * @param E_kin - variable to save the kinetic energy
 */
void VMC::E_Pot_Kin(double alpha, double beta, double& E_pot, double& E_kin) {

    E_pot = H->H_potential(Trial);
    E_kin = H->H_kinetic(Trial, alpha, beta);
}

/**
 * Update of the statistics for each Monte Carlo sample
 * @param del_E  - contribution of local energy
 */
void VMC::update_statistics(double del_E) {

    E += del_E;
    E_sq += del_E*del_E;

#if MINIMIZE

    exp_par_psi += par_psi;
    exp_par_psi2 += del_E*par_psi;

#elif DISTANCE

    for (int j = 1; j < numpart; j++) {
        for (int i = 0; i < j; i++) {
            r_dist += Trial->Pos->r_int(i, j);
            r_distsq += Trial->Pos->r_int(i, j) * Trial->Pos->r_int(i, j);
        }
    }

#endif  

    return;

}

/**
 * Update of the statistics if E_POT_KIN = "true"
 * @param del_E -  total local energy
 * @param del_Epot - contribution of potential energy
 * @param del_Ekin - contribution of kinetic energy
 */
void VMC::update_statistics(double del_E, double del_Epot, double del_Ekin) {
    E += del_E;
    E_sq += del_E*del_E;
    E_pot += del_Epot;
    E_potsq += del_Epot*del_Epot;
    E_kin += del_Ekin;
    E_kinsq += del_Ekin*del_Ekin;

#if MINIMIZE

    exp_par_psi += par_psi;
    exp_par_psi2 += del_E*par_psi;

#elif DISTANCE

    for (int j = 1; j < numpart; j++) {
        for (int i = 0; i < j; i++) {
            r_dist += Trial->Pos->r_int(i, j);
            r_distsq += Trial->Pos->r_int(i, j) * Trial->Pos->r_int(i, j);
        }
    }

#endif  

    return;

}

/**
 * If the move is not accepted, reset position and Slater inverse
 * @param p - particle that has been moved
 */
void VMC::not_accept(int p) {

    for (int k = 0; k < dim; k++) {
        R_tr(p, k) = R_cur(p, k);
    }

    Trial->Pos_tr->update_p(p, R_cur);
    Trial->SlaterPsi->reset(p);

    return;

}

/**
 * Compute numerically the derivative of the wavefunction with respect to
 * the variational parameters
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void VMC::part_psi(double alpha, double beta) {

    double h = 0.002;
    mat R = Trial->Pos->current;

    double al_plus = Trial->value(R, alpha + h, beta);
    double al_minus = Trial->value(R, alpha - h, beta);
    double beta_plus = Trial->value(R, alpha, beta + h);
    double beta_minus = Trial->value(R, alpha, beta - h);

    // Resets also internal values in objects
    double wf_old = Trial->value(R, alpha, beta);

    double deriv_alpha = (al_plus - al_minus) / (2 * h);
    double deriv_beta = (beta_plus - beta_minus) / (2 * h);

    par_psi(0) = deriv_alpha;
    par_psi(1) = deriv_beta;

    par_psi /= wf_old;

    return;

}

/**
 * Compute the derivative of the wavefunction with respect to alpha 
 * numerically and with respect to beta analytically
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void VMC::part_psi_beta_analytic(double alpha, double beta) {

    double rij;
    par_psi.zeros();
    
    // analytical derivative with respect to beta
    for(int j = 0; j<numpart; j++){
        for(int i = 0; i<j; i++){
            rij = Trial->Pos->r_int(i,j);
            par_psi(1) -= Trial->JastrowPsi->a(i,j)*rij*rij/
                    ((1+beta*rij)*(1+beta*rij));
        }
    }
    
    // numerical derivative with respect to beta
    double h = 0.005;
    mat R = Trial->Pos->current;

    double al_plus = Trial->value(R, alpha + h, beta);
    double al_minus = Trial->value(R, alpha - h, beta);
    
    // Resets also internal values in objects
    double wf_old = Trial->value(R, alpha, beta);

    double deriv_alpha = (al_plus - al_minus) / (2 * h);
    
    par_psi(0) = deriv_alpha;
    par_psi /= wf_old;
    
    return;
    
}

/**
 * Compute average distance between two particles
 * @return distance
 */
double VMC::get_rdist() {

    int num_intact; // number of interactions
    num_intact = (numpart * numpart - numpart) / 2;

    r_dist /= num_intact;
    return r_dist;
}

/*
 * Returns  the average of the squared distance between two particles, 
 * Needed to compute the statistical error of <r>
 */
double VMC::get_rdistsq() {

    int num_intact; // number of interactions
    num_intact = (numpart * numpart - numpart) / 2;

    r_distsq /= num_intact;
    return r_distsq;
}

/**
 * @return energy expectation value
 */
double VMC::get_E() const {
    return E;
}

/** 
 * @return expectation value of the square of the energy
 */
double VMC::get_Esq() const {
    return E_sq;
}

/**
 * @return expectation value of the potential energy
 */
double VMC::get_Epot() const {
    return E_pot;
}

/**
 * @return expectation value of the square of the potential energy
 */
double VMC::get_Epotsq() const {
    return E_potsq;
}

/**
 * @return expectation value of the kinetic energy
 */
double VMC::get_Ekin() const {
    return E_kin;
}

/**
 * @return expectation value of the square of the kinetic energy
 */
double VMC::get_Ekinsq() const {
    return E_kinsq;
}

/**
 * @return first vector needed for CGM
 */
vec VMC::get_par_psi() const {
    return exp_par_psi;
}

/**
 * @return second vector needed for CGM
 */
vec VMC::get_par_psi2() const {
    return exp_par_psi2;
}

/**
 * The VMC algorithm adapted to the blocking procedure
 * @param N - Number of MC cycles
 * @param N_therm - Number of thermalization steps
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return local energies of all samples
 */
vec VMC::run_algo_blocking(int N, int N_therm, double alpha, double beta) {

    int i;
    double rat, eps;
    double del_E = 0;
    vec all_energies(N * numpart);


    // ************************* Thermalization ******************************

    accepted = 0;

    initialize(alpha, beta); // Initialize the system

    for (i = 0; i < N_therm; i++) {

        for (int p = 0; p < numpart; p++) {// Loop over all particles

            // Calculate trial position
            trial_pos(p, alpha, beta);

            // Compute acceptance ratio
            rat = ratio(p, alpha, beta);

            // Check if move is accepted
            if (rat >= 1.0) { // accept if probability is greater
                accept(p, alpha, beta);
                accepted++;
            } else { // otherwise check against random number

                eps = ran3(&idum);
                if (eps < rat) {
                    accept(p, alpha, beta);
                    accepted++;
                } else
                    not_accept(p); // Do not accept
            }
        }
    }


    // *********************** After thermalization  **************************

    if (N > 0) accepted = 0; // For function delta_opt()

    for (i = 0; i < N; i++) {

        for (int p = 0; p < numpart; p++) { // Loop over all particles

            // Calculate trial position
            trial_pos(p, alpha, beta);

            // Compute acceptance ratio
            rat = ratio(p, alpha, beta);

            // Check if move is accepted
            if (rat >= 1.0) { // accept if probability is greater
                accept(p, alpha, beta);
                del_E = E_local(alpha, beta);
                accepted++;
            } else { // otherwise check against random number

                eps = ran3(&idum);
                if (eps < rat) {
                    accept(p, alpha, beta);
                    del_E = E_local(alpha, beta);
                    accepted++;
                } else
                    not_accept(p);

            }

            // Collecting the local energies for blocking
            all_energies(i) = del_E;

        }
    }

    return all_energies;

}



////////////////////////////////////////////////////////////////////////////////
//                           class Metropolis

/** @brief Class for Metropolis algorithm, subclass of VMC
    @author sarahrei
    @date 11 April  2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param Trial - pointer to trial wavefunction
 * @param H - pointer to Hamiltonian for the local energy
 * @param idum - seed for random number generator
 */
Metropolis::Metropolis(Wavefunction* Trial, Hamiltonian* H, long idum) {

    this->Trial = Trial;
    this->dim = Trial->get_dim();
    this->numpart = Trial->get_numpart();

    this->H = H;
    this->idum = idum;

    R_cur.set_size(numpart, dim);
    R_tr.set_size(numpart, dim);
    exp_par_psi = zeros(2);
    exp_par_psi2 = zeros(2);
    par_psi = zeros(2);

    interaction = Trial->get_int();

}

/**
 * Computation of the trial position according to the Metropolis algorithm
 * @param p - particle that has been moved compared to previous configuration
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void Metropolis::trial_pos(int p, double alpha, double beta) {

    for (int k = 0; k < dim; k++) {
        R_tr(p, k) = R_cur(p, k) + delta * (2 * ran3(&idum) - 1);
    }

    Trial->Pos_tr->update_p(p, R_tr);

    return;
}

/**
 * Compute acceptance ratio according to the Metropolis algorithm
 * @param p - particle that has been moved 
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return 
 */
double Metropolis::ratio(int p, double alpha, double beta) {

    double r;

    r = Trial->ratio(R_tr, p, alpha, beta);
    r *= r;

    return r;
}

/**
 * Initializations for the Metropolis algorithm
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void Metropolis::initialize(double alpha, double beta) {

    E = 0.0;
    E_sq = 0.0;
    E_kin = 0.0;
    E_kinsq = 0.0;
    E_pot = 0.0;
    E_potsq = 0.0;
    exp_par_psi.zeros(),
            exp_par_psi2.zeros();
    r_dist = 0.0;
    r_distsq = 0.0;

    // Initialize the position
    R_cur.randu();
    R_tr = R_cur;
    Trial->Pos->update(R_cur);
    Trial->Pos_tr->update(R_cur);

    // Initialize Slater determinant and Jastrow factor
    Trial->SlaterPsi->value(R_cur, alpha);
    if (interaction) Trial->JastrowPsi->value(Trial->Pos, beta);

    // Initialize Slater inverse
    int n2 = numpart / 2;
    for (int i = 0; i < 2; i++) { // For spin up and down
        Trial->SlaterPsi->slat_inv.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2)
                = inv(Trial->SlaterPsi->slater.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2));
    }

    return;
}

/**
 * Performs all updates if a move has been accepted
 * @param p - particle that has been moved 
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void Metropolis::accept(int p, double alpha, double beta) {


    for (int k = 0; k < dim; k++) {
        R_cur(p, k) = R_tr(p, k); // Update position
    }

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
 * Function to determine the step length delta. Aim is to accept about 
 * 50% of the moves. The algorithm is based on bisection, i.e. start 
 * values for minimum and maximum have to be provided, as well as the 
 * tolerance interval "eps". 
 * @param N_delta - number of MC steps per trial of delta
 * @param alpha - first variational parameter 
 * @param beta - second variational parameter
 * @param del_min - minimal value of delta
 * @param del_max - maximal value of delta
 * @param eps - tolerance interval between del_min and del_max
 */
void Metropolis::delta_opt(int N_delta, double alpha, double beta, double del_min,
        double del_max, double eps) {

    double accept_min; // rate of acceptance for del_min
    double accept_mid; // rate of acceptance for mid point              

    while ((del_max - del_min) > eps) {

        delta = del_min;
        run_algo(0, N_delta, alpha, beta);
        accept_min = (accepted * 1.0) / (N_delta * numpart);

        delta = (del_max + del_min) / 2.;
        run_algo(0, N_delta, alpha, beta);

        accept_mid = (accepted * 1.0) / (N_delta * numpart);

        // Update deLmin and del_max
        if ((accept_min - 0.5)*(accept_mid - 0.5) < 0) {
            del_max = (del_max + del_min) / 2;
        } else {
            del_min = delta;
        }

    } // End while loop

    delta = (del_max + del_min) / 2.;

    return;
}



////////////////////////////////////////////////////////////////////////////////
//                      class Metropolis-Hastings

/** @brief Class for Metropolis-Hastings algorithm, subclass of VMC
    @author sarahrei
    @date 11 April  2012
 */
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param Trial - pointer to trial wavefunction
 * @param H - pointer to Hamiltonian for the local energy
 * @param seed - seed for normal random number generator
 * @param idum - seed for uniform random number generator
 */
Metropolis_Hastings::Metropolis_Hastings(Wavefunction* Trial, Hamiltonian* H,
        long idum) {

    this->Trial = Trial;
    this->dim = Trial->get_dim();
    this->numpart = Trial->get_numpart();
    this->H = H;
    this->idum = idum;
    
    QFo = new QForce(Trial);

    R_cur.set_size(numpart, dim);
    R_tr.set_size(numpart, dim);
    delta = 0.5;

    interaction = Trial->get_int();

    seed = (int) abs(idum);
    RanNormalSetSeedZigVec(&seed, 200); // Initialize random number generator

    exp_par_psi = zeros(2);
    exp_par_psi2 = zeros(2);
    par_psi = zeros(2);

}

/**
 * Compute acceptance ratio according to the Metropolis-Hastings algorithm
 * @param p - particle that has been moved 
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return 
 */
void Metropolis_Hastings::trial_pos(int p, double alpha, double beta) {

    double chi;
    bruteforce_tour = false;

    // Compute trial position Metropolis-Hastings
    for (int j = 0; j < dim; j++) {
        chi = DRanNormalZigVec();
        R_tr(p, j) = R_cur(p, j) + 0.5 * QFo->qf_old(p, j) * del_t + chi * sqrt(del_t);

    }

    Trial->Pos_tr->update_p(p, R_tr);
    test_inv = Trial->SlaterPsi->inverse(p, R_tr, alpha);

    // New quantum force
    QFo->new_qf(p, R_tr, Trial->Pos_tr, alpha, beta, test_inv);

    // In case of divergencies, take brute force step
    if (QFo->qf_new.max() > 100 || QFo->qf_new.min() < -100)
        bruteforce_tour = true;

    if (bruteforce_tour) {

        for (int j = 0; j < dim; j++) {
            chi = DRanNormalZigVec();
            R_tr(p, j) = R_cur(p, j) + delta * (2 * ran3(&idum) - 1);
        }

        Trial->Pos_tr->update_p(p, R_tr);
        test_inv = Trial->SlaterPsi->inverse(p, R_tr, alpha);
        QFo->new_qf(p, R_tr, Trial->Pos_tr, alpha, beta, test_inv);
    }

    return;

}

/**
 * Compute acceptance ratio according to the Metropolis-Hastings algorithm
 * @param p - particle that has been moved 
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return 
 */
double Metropolis_Hastings::ratio(int p, double alpha, double beta) {

    double r, ef;
    vec ef1(dim), ef2(dim);

    r = Trial->ratio(R_tr, p, alpha, beta);
    r *= r;

    if (bruteforce_tour) return r;

    for (int k = 0; k < dim; k++) {
        ef1(k) = QFo->qf_old(p, k) + QFo->qf_new(p, k);
        ef2(k) = (QFo->qf_old(p, k) - QFo->qf_new(p, k))*0.25 * del_t;
        ef2(k) += R_cur(p, k) - R_tr(p, k);
    }

    ef = 0.5 * dot(ef1, ef2);
    ef = exp(ef);
    r *= ef;

    return r;

}

/**
 * Initializations for the Metropolis-Hastings algorithm
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void Metropolis_Hastings::initialize(double alpha, double beta) {

    E = 0.0;
    E_sq = 0.0;
    E_kin = 0.0;
    E_kinsq = 0.0;
    E_pot = 0.0;
    E_potsq = 0.0;
    exp_par_psi.zeros(),
            exp_par_psi2.zeros();
    r_dist = 0.0;
    r_distsq = 0.0;

    // Initialize the position
    R_cur.randu();
    R_tr = R_cur;
    Trial->Pos->update(R_cur);
    Trial->Pos_tr->update(R_cur);

    // Initialize Slater determinant and Jastrow factor
    Trial->SlaterPsi->value(R_cur, alpha);
    Trial->JastrowPsi->value(Trial->Pos, beta);

    // Initialize Slater inverse
    int n2 = numpart / 2;

    for (int i = 0; i < 2; i++) { // For spin up and down
        Trial->SlaterPsi->slat_inv.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2)
                = inv(Trial->SlaterPsi->slater.submat(0, i*n2, n2 - 1, n2 - 1 + i * n2));
    }
    test_inv = Trial->SlaterPsi->slat_inv;

    for (int p = 0; p < numpart; p++) {
        QFo->new_qf(p, R_cur, Trial->Pos, alpha, beta, test_inv);
    }

    QFo->qf_old = QFo->qf_new;

    return;
}

/**
 * Performs all updates if a move has been accepted
 * @param p - particle that has been moved 
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 */
void Metropolis_Hastings::accept(int p, double alpha, double beta) {


    for (int k = 0; k < dim; k++) {
        R_cur(p, k) = R_tr(p, k); // Update position
        QFo->qf_old(p, k) = QFo->qf_new(p, k);
    }


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
 * Set timestep delta_t
 * @param del_t - timestep  delta_t
 */
void Metropolis_Hastings::set_delt(double del_t) {
    this->del_t = del_t;

}


