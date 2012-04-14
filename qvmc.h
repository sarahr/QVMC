/* 
 * File:   qvmc.h
 * Author: sarahrei 
 *
 * Created on 19. januar 2012, 12:14 
 */

#ifndef QVMC_H
#define	QVMC_H 

#include "Radial.h"
#include "Hermite.h"
#include "Jastrow.h"
#include "Slater.h"
#include "ExpFactor.h"
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "QForce.h"

//#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>
#include <armadillo> 

using namespace std;
using namespace arma;

/** @brief Class for the VMC algorithm
    @author sarahrei
    @date 11 April  2012
 */
class VMC { 
protected:
    int numpart;
    int dim;
    Wavefunction* Trial;
    Hamiltonian* H;

    double delta;
    double E, E_sq, E_kin, E_kinsq, E_pot, E_potsq;
    vec exp_par_psi, exp_par_psi2, par_psi;
    double accepted;
    double r_dist, r_distsq;

    long idum;
    mat R_cur, R_tr;
    bool interaction;

public:

    VMC() {
    };

    /**
     * Constructor
     * @param N - number of MC cycles
     * @param N_therm - number of thermalization steps
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void run_algo(int N, int N_therm, double alpha, double beta);

    /**
     * Compute the local energy
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return local energy
     */
    double E_local(double alpha, double beta);

    /**
     * Compute kinetic and potential energy separately
     * @param alpha - first variational parameter
     * @param beta  - second variational parameter
     * @param E_pot - variable to save the potential energy
     * @param E_kin - variable to save the kinetic energy
     */
    void E_Pot_Kin(double alpha, double beta, double& E_pot, double& E_kin);

    /**
     * Update of the statistics for each Monte Carlo sample
     * @param del_E  - contribution of local energy
     */
    void update_statistics(double del_E);
    
    /**
     * Update of the statistics if E_POT_KIN = "true"
     * @param del_E -  total local energy
     * @param del_Epot - contribution of potential energy
     * @param del_Ekin - contribution of kinetic energy
     */
    void update_statistics(double del_E, double del_Epot, double del_Ekin);

    /**
     * Compute numerically the derivative of the wavefunction with respect to
     * the variational parameters
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void part_psi(double alpha, double beta);

    /**
     * If the move is not accepted, reset position and Slater inverse
     * @param p - particle that has been moved
     */
    void not_accept(int p);

    /**
     * The VMC algorithm adapted to the blocking procedure
     * @param N - Number of MC cycles
     * @param N_therm - Number of thermalization steps
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return local energies of all samples
     */
    vec run_algo_blocking(int N, int N_therm, double alpha, double beta);

    /**
     * Virtual function for trial position
     * @param p - parameter that has been moved
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    virtual void trial_pos(int p, double alpha, double beta) = 0;

    /**
     * virtual function for acceptance ratio
     * @param p - parameter that has been moved
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return acceptance ratio
     */
    virtual double ratio(int p, double alpha, double beta) = 0;

    /**
     * virtual function for initializations of VMC algorithm
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    virtual void initialize(double alpha, double beta) = 0;

    /**
     * Virtual function for updating after accepted Metropolis steps
     * @param p - parameter that has been moved
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    virtual void accept(int p, double alpha, double beta) = 0;

    /**
     * @return energy expectation value
     */
    double get_E() const;

    /**
     * @return expectation value of the square of the energy
     */
    double get_Esq() const;

    /**
     * @return expectation value of potential energy
     */
    double get_Epot() const;

    /**
     * @return expectation value of the square of the potential energy
     */
    double get_Epotsq() const;

    /**
     * @return expectation value of kinetic energy
     */
    double get_Ekin() const;

    /**
     * @return expectation value of the square of the kinetic energy
     */
    double get_Ekinsq() const;

    /**
     * @return first vector needed for CGM
     */
    vec get_par_psi() const;

    /**
     * @return second vector needed for CGM
     */
    vec get_par_psi2() const;

    /**
     * Compute average distance between two particles
     * @return average distance
     */
    double get_rdist();

    /**
     * Compute  the average of the squared distance between two particles, 
     * Needed to compute the statistical error of <r>
     * @return average of the squared distance
     */
    double get_rdistsq();
};

/** @brief Class for Metropolis algorithm, subclass of VMC
    @author sarahrei
    @date 11 April  2012
 */
class Metropolis : public VMC {
public:

    /**
     * Constructor
     * @param Trial - pointer to trial wavefunction
     * @param H - pointer to Hamiltonian for the local energy
     * @param idum - seed for random number generator
     */
    Metropolis(Wavefunction* Trial, Hamiltonian* H, long idum);

    /**
     * Computation of the trial position according to the Metropolis algorithm
     * @param p - particle that has been moved compared to previous configuration
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void trial_pos(int p, double alpha, double beta);

    /**
     * Compute acceptance ratio according to the Metropolis algorithm
     * @param p - particle that has been moved 
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return acceptance ratio according to the Metropolis algorithm
     */
    double ratio(int p, double alpha, double beta);

    /**
     * Initializations for the Metropolis algorithm
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void initialize(double alpha, double beta);

    /**
     * Performs all updates if a move has been accepted
     * @param p - particle that has been moved 
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void accept(int p, double alpha, double beta);

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
    void delta_opt(int N_delta, double alpha, double beta, double del_min,
            double del_max, double eps);

};

/** @brief Class for Metropolis-Hastings algorithm, subclass of VMC
    @author sarahrei
    @date 11 April  2012
 */
class Metropolis_Hastings : public VMC {
private:
    QForce* QFo;
    double del_t;
    int seed;
    mat test_inv;
    int count0;
    bool bruteforce_tour;

public:

    /**
     * Constructor
     * @param Trial - pointer to trial wavefunction
     * @param H - pointer to Hamiltonian for the local energy
     * @param seed - seed for normal random number generator
     * @param QFo - pointer to quantum force
     * @param idum - seed for uniform random number generator
     */
    Metropolis_Hastings(Wavefunction* Trial, Hamiltonian* H, int seed, QForce* QFo, long idum);

    /**
     * Computation of the trial position according to the Metropolis-Hastings algorithm
     * @param p - particle that has been moved compared to previous configuration
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void trial_pos(int p, double alpha, double beta);

    /**
     * Compute acceptance ratio according to the Metropolis-Hastings algorithm
     * @param p - particle that has been moved 
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return acceptance ratio according to the Metropolis-Hastings algorithm
     */
    double ratio(int p, double alpha, double beta);

    /**
     * Initializations for the Metropolis-Hastings algorithm
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void initialize(double alpha, double beta);

    /**
     * Performs all updates if a move has been accepted
     * @param p - particle that has been moved 
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     */
    void accept(int p, double alpha, double beta);

    /**
     * Set timestep delta_t
     * @param del_t - timestep delta_t
     */
    void set_delt(double del_t);

};

#endif	/* QVMC_H */
