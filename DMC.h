/* 
 * File:   DMC.h
 * Author: sarahrei
 *
 * Created on April 22, 2012, 7:13 PM
 */

#ifndef DMC_H
#define	DMC_H

#include "Walker.h"
#include "ControlWalkers.h"

#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>
//#include <armadillo> 

using namespace std;
using namespace arma;

/** @brief Class for the DMC algorithm
    @author sarahrei
    @date 15 May 2012
 */
class DMC {
private:

    int numpart;
    int dim;
    double alpha, beta;

    Walker** QDot;
    ControlWalkers* Authority;
    Hamiltonian* H;

    long idum;
    double E_old, E_new;
    int jastrow;
    double omega;
    double del_t;

public:


    /**
     * Constructor
     * @param numpart - number of particles
     * @param dim - dimension
     * @param omega - oscillator frequency
     * @param jastrow  - 0: without, 1: with Jastrow factor
     * @param H - Hamiltonian of the system
     * @param alpha  - first variational parameter of trial wavefunction
     * @param beta - second variational parameter of trial wavefunction
     * @param myrank - rank in MPI
     */
    DMC(int numpart, int dim, double omega, int jastrow, Hamiltonian* H,
            double alpha, double beta, int myrank);

    /**
     * Destructor
     */
    ~DMC();


    /**
     * Main DMC algorithm
     * @param N - number of cycles per block in sampling phase
     * @param N_therm - number of thermalization steps
     * @param del_t- timestep
     * @param E_t - reference energy from VMC
     * @param samples - number of samples/blocks in sampling phase
     * @param N_equi - number of cycles per block in equilibration phase
     * @param equi_samples - number of samples for equilibration
     * @param nw_init - initial/target number of walkers
     * @param outfile - name of outputfile
     * @param myrank - rank in MPI
     * @param numprocs - number of MPI processes
     */
    void run_algo(int N, int N_therm, double del_t, double E_t, int samples,
            int N_equi, int equi_samples, int nw_init, char* outfile, int myrank, int numprocs);
     

    /**
     * Thermalization period to get a representative start position for all walkers
     * Movement according to Metropolis-Hastings algorithm
     * @param N_therm - number thermalization steps
     * @param e_local_old - array with local energies
     * @param first - first walker to be moved
     * @param last - (last-1) walker to be moved
     */
    void initialize(int N_therm, double *e_local_old, int nw_init);

    /**
     * Sorting of the walkers such that the array of walkers is dense again,
     * "num_alive" is updated
     * @param num_alive - number of walkers that are alive
     * @param killed - number of walkers that have been killed
     * @param num_copies - number of new walker-copies
     * @param alive - array storing alive/dead for each walker
     * @param e_local_old - array with local energies from previous cycle
     */
    void sortWalkers(int &num_alive, int killsd, int num_resurrected, bool *alive, double *e_local_old);


};


#endif	/* DMC_H */

