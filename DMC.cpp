
#include "DMC.h" 
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <time.h>
#include "lib.h"
#include <mpi.h>


#define BLOCKING false
#define DENSITY true // store equilibrium walker distribution


////////////////////////////////////////////////////////////////////////////////
//                             class DMC

/** @brief Class for the DMC algorithm
    @author sarahrei
    @date 15 May 2012
 */
////////////////////////////////////////////////////////////////////////////////

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
DMC::DMC(int numpart, int dim, double omega, int jastrow, Hamiltonian* H,
        double alpha, double beta, int myrank) {

    this->numpart = numpart;
    this->dim = dim;
    this->omega = omega;
    this->jastrow = jastrow;
    this->alpha = alpha;
    this->beta = beta;
    this->H = H;

    Authority = new ControlWalkers(numpart, dim);

    // Initialize random number generators
    int seed = time(NULL) + myrank * 10;
    int seed2 = seed + myrank * 3;
    RanNormalSetSeedZigVec(&seed, 200);
    RanSetSeed_MWC8222(&seed2, 100);

}

/**
 * Destructor
 */
DMC::~DMC() {

    delete [] QDot;
    delete Authority;
}

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
void DMC::run_algo(int N, int N_therm, double del_t, double E_t, int samples,
        int N_equi, int equi_samples, int nw_init, char* outfile, int myrank, int numprocs) {

    ofstream ofile;
    ofile.open(outfile);

#if BLOCKING
    ofstream ofile2;
    ostringstream ost;
    ost << "blocks_rank" << myrank << ".dat";
    ofile2.open(ost.str().c_str(), ios::out | ios::binary);
#endif
    
#if DENSITY
    ofstream ofile3;
    ostringstream ost;
    ost << "density" << myrank << ".dat";
    ofile3.open(ost.str().c_str(), ios::out);
#endif

    double e_local = 0.0;
    double e_local_temp = 0.0;
    int i, j;
    int num_alive = nw_init;
    double green_branch;
    int num_next_gen;
    int num_copies = 0;
    int killed = 0;
    double rat, eps, elapsed_time;
    bool mh_test;
    double E_old, E_new, E_new_total;
    int sorted_in_this_block, sorted_total, equi;
    int total_walkers;
    double e_ref = E_t;

    // Create array of walkers    
    QDot = new Walker*[nw_init * 3];
    bool* alive = new bool[nw_init * 3];
    double* e_local_old = new double[nw_init * 3];

    // Initialize the walkers  
    for (i = 0; i < nw_init * 3; i++) {

        QDot[i] = new Walker(numpart, dim, omega, jastrow, H, alpha, beta);
        QDot[i]->initWalker(del_t);

        if (i < nw_init) alive[i] = true;
        else alive[i] = false;

    }

    // Thermalization of walkers according to Metropolis algorithm
    initialize(N_therm, e_local_old, nw_init);

    if (myrank == 0) cout << "Thermalization finished.\n";
    fflush(stdout);

    E_old = E_t;
    E_new = E_t;
    num_alive = nw_init;

    ///////////////////////////////////////////////////////////////////////////
    //                        Equilibration phase
    ///////////////////////////////////////////////////////////////////////////
    
    equi = 0;

    while (equi < equi_samples) {

        equi++;
        E_old = E_new;
        total_walkers = 0;
        e_local = 0.0;
        e_local_temp = E_t;
        sorted_in_this_block = 0;

        for (int k = 0; k < N_equi; k++) {

            num_copies = 0;
            killed = 0;
            int all_counts = 0;

            for (int loop_p = 0; loop_p < num_alive; loop_p++) {

                // Move all particles, one at a time
                for (int p = 0; p < numpart; p++) {

                    // Calculate trial position
                    QDot[loop_p]->trial_pos(p);

                    // Compute acceptance ratio
                    rat = QDot[loop_p]->ratio(p);

                    // Check if move is accepted
                    eps = DRan_MWC8222();

                    if (eps < rat)
                        mh_test = true;
                    else mh_test = false;

                    // Accept if Metropolis test accepted and no node crossed
                    if (mh_test && !QDot[loop_p]->nodeCrossed()) {

                        QDot[loop_p]->accept(p);

                    } else
                        QDot[loop_p]->not_accept(p);
                } 

                // Compute local energy
                 e_local_temp = QDot[loop_p]->E_local();

                // Calculate new branching factor
                green_branch = exp(-del_t * (0.5 * (e_local_temp + 
                        e_local_old[loop_p]) - e_ref));

                // Compute number of copies in next generation
                num_next_gen = (int) (green_branch + DRan_MWC8222());
                
                // Prevent explosions
                if (num_next_gen > 10) num_next_gen = 10;

                // Energy cutoff to avoid divergencies near nodes
                if (e_local_temp < E_t - 2 / sqrt(del_t)) {
                    e_local_temp = E_t - 2 / sqrt(del_t);
                    num_next_gen = 0;
                } else if (e_local_temp > E_t + 2 / sqrt(del_t)) {
                    e_local_temp = E_t + 2 / sqrt(del_t);
                    num_next_gen = 0;
                }

                
                if(green_branch > 10)  green_branch = 10;
                
                // Update local energy weighted with branching factor
                e_local += green_branch*e_local_temp;
                all_counts++;

                // Branching
                if (num_next_gen == 0) {
                    alive[loop_p] = false;
                    killed++;

                } else if (num_next_gen >= 2) {

                    j = num_alive + num_copies;

                    while (num_next_gen > 1) {
                        if (j < 3 * nw_init) {
                            Authority->copyWalker(QDot[loop_p], QDot[j]);
                            e_local_old[j] = e_local_temp;
                            alive[j] = true;
                            num_copies++;
                            j++;
                        }
                        num_next_gen--;
                    }
                }
                 e_local_old[loop_p] = e_local_temp;

            }


            total_walkers += all_counts;
            sorted_in_this_block = sorted_in_this_block - killed + num_copies;
            
            // Sort the walkers
            sortWalkers(num_alive, killed, num_copies, alive, e_local_old);

            // Update reference energy
            e_ref = E_old - log((double) num_alive / (double) nw_init) ;
            
        }
      
        e_local /= total_walkers;
        E_new = e_local;
        cout << num_alive << endl;

        // Data exchange between MPI processes
        MPI_Reduce(&sorted_in_this_block, &sorted_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&E_new, &E_new_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (myrank == 0) {
            E_new_total /= numprocs;
            cout << "Sorted: " << sorted_total << endl;
            cout << setprecision(8) << E_new_total << endl;
            ofile << equi * N_equi << " " << sorted_in_this_block << E_new << endl;
        }

    }
    
    //////////////////////////////////////////////////////////////////////////
    //                        Sampling phase
    //////////////////////////////////////////////////////////////////////////


    if (myrank == 0) {
        cout << "Begin sampling\n";
        ofile << endl;
    }

    double ener_cumu = 0.0;
    int no_sample = 0;
    double e_floating = 0.0;
    vec E_i(samples);
    int dist = samples*N/10000;

    while (no_sample < samples) {

        E_old = E_new;
        int total_walkers = 0;
        e_local = 0.0;
        sorted_in_this_block = 0;
        double sum_sq = 0.0;

        for (int k = 0; k < N; k++) {

            num_copies = 0;
            killed = 0;
            int all_counts = 0;
            
            // Loop over walkers
            for (int loop_p = 0; loop_p < num_alive; loop_p++) {

                // Move all particles, one at a time
                for (int p = 0; p < numpart; p++) {

                    // Calculate trial position
                    QDot[loop_p]->trial_pos(p);

                    // Compute acceptance ratio
                    rat = QDot[loop_p]->ratio(p);

                    // Check if move is accepted
                    eps = DRan_MWC8222();

                    if (eps < rat)
                        mh_test = true;
                    else mh_test = false;

                    // Reject walkers that have crossed node or fail the MH_test
                    if (mh_test && !QDot[loop_p]->nodeCrossed()) {
                        QDot[loop_p]->accept(p);
  
                    } else
                        QDot[loop_p]->not_accept(p);
                }

                //  Find local energy
                e_local_temp = QDot[loop_p]->E_local();

                // Compute new branching factor
                green_branch = exp(-del_t * (0.5 * (e_local_temp + e_local_old[loop_p]) - e_ref));
                
                // Compute number of copies in next generation
                num_next_gen = (int) (green_branch + DRan_MWC8222());

                // Prevent explosions
                if (num_next_gen > 10) num_next_gen = 10;

                // Energy cutoff
                if (e_local_temp < E_t - 2 / sqrt(del_t)) {
                    e_local_temp = E_t - 2 / sqrt(del_t);
                    num_next_gen = 0;
                } else if (e_local_temp > E_t + 2 / sqrt(del_t)) {
                    e_local_temp = E_t + 2 / sqrt(del_t);
                    num_next_gen = 0;
                }
                     
                
                if(green_branch > 10)  green_branch = 10;
                if(green_branch*e_local_temp > 100000 || green_branch*e_local_temp < -100000 ) cout << "Divergenz" << endl;
                e_local += green_branch*e_local_temp;
                sum_sq += e_local_temp * e_local_temp * green_branch*green_branch;
                all_counts++;
#if BLOCKING
                ofile2 << e_local_temp << endl;
#endif
                
#if DENSITY
                if (k % dist == 0) {
                    for (int l = 0; l < numpart; l++) {
                        ofile3 << sqrt(QDot[loop_p]->Trial->Pos->r(l)) << " ";
                    }
                    
                    ofile3 << green_branch << endl;
                }
#endif
                
                // Branching
                if (num_next_gen == 0) {
                    alive[loop_p] = false;
                    killed++;

                } else if (num_next_gen >= 2) {

                    j = num_alive + num_copies;

                    while (num_next_gen > 1) {
                        if (j < 3 * nw_init) {
                            Authority->copyWalker(QDot[loop_p], QDot[j]);
                            e_local_old[j] = e_local_temp;
                            alive[j] = true;
                            num_copies++;
                            j++;
                        }
                        num_next_gen--;
                    }
                }
                e_local_old[loop_p] = e_local_temp;

            }

            total_walkers += all_counts;
            sorted_in_this_block = sorted_in_this_block - killed + num_copies;
            
            // Sort the walkers
            sortWalkers(num_alive, killed, num_copies, alive, e_local_old);

        }

        e_local /= total_walkers;
        E_new = e_local;
        e_ref = E_old;

        // Communication between MPI processes
        MPI_Reduce(&sorted_in_this_block, &sorted_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&E_new, &E_new_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (myrank == 0) {

            E_new_total /= numprocs;
            cout << "Sorted: " << sorted_total << endl;

            E_i(no_sample) = E_new_total;
            ener_cumu += E_new_total;
            no_sample++;
            e_floating = ener_cumu / no_sample;

            cout << E_new_total << " " << e_floating << endl;
            ofile << equi * N_equi + no_sample * N << " " << sorted_total << " "
                    << setprecision(12) << E_new_total << " " << e_floating << endl;

        } else
            no_sample++;

    } // End all samples

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0) {
        
        double var2 = 0;
        for (int k = 0; k < samples; k++) {
            var2 += E_i(k) * E_i(k);
        }

        var2 /= samples;
        var2 -= ener_cumu / samples * ener_cumu / samples;

        cout << "\nFinal energy: " << setprecision(16) << ener_cumu / samples << "\n";
        ofile << "\nFinal energy: " << setprecision(16) << ener_cumu / samples << "\n";
        cout << "\nVariance: " << setprecision(16) << var2 << " " << sqrt(var2 / samples) << "\n";
        ofile << "\nVariance: " << setprecision(16) << var2 << " " << sqrt(var2 / samples) << "\n";
        ofile.close();
    }

     MPI_Barrier(MPI_COMM_WORLD);

#if BLOCKING
    ofile2.close();
#endif   
#if DENSITY
    ofile3.close();
#endif

}


/**
 * Thermalization period to get a representative start position for all walkers
 * Movement according to Metropolis-Hastings algorithm
 * @param N_therm - number thermalization steps
 * @param e_local_old - array with local energies
 * @param first - first walker to be moved
 * @param last - (last-1) walker to be moved
 */
void DMC::initialize(int N_therm, double *e_local_old, int nw_init) {

    double rat, eps;
    int i, j;
    double del_E;

    // Loop over walkers
    for (i = 0; i < nw_init; i++) {

        // Loop over themalization steps
        for (j = 0; j < N_therm; j++) {

            //  Move all particles, one at a time
            for (int p = 0; p < numpart; p++) {

                // Calculate trial position
                QDot[i]->trial_pos(p);

                // Compute acceptance ratio
                rat = QDot[i]->ratio(p);

                // Check if move is accepted
                eps = DRan_MWC8222();

                if (eps < rat) {
                    QDot[i]->accept(p); // Accept the move

                } else {
                    QDot[i]->not_accept(p); // Do not accept
                }


            } // All particles moved

            if (j == N_therm - 1) {
                e_local_old[i] = QDot[i]->E_local();
            }
        }
    }

}

/**
 * Sorting of the walkers such that the array of walkers is dense again,
 * "num_alive" is updated
 * @param num_alive - number of walkers that are alive
 * @param killed - number of walkers that have been killed
 * @param num_copies - number of new walker-copies
 * @param alive - array storing alive/dead for each walker
 * @param e_local_old - array with local energies from previous cycle
 */
void DMC::sortWalkers(int &num_alive, int killed, int num_copies, bool *alive,
        double *e_local_old) {

    int sorted = 0, i = 0;
    int j = num_alive + num_copies - 1; // Starting at the end of list

    while (sorted < killed) { // Fill all places where walkers died
        
        Walker* temp;

        if (alive[i] == false) {

            while (!alive[j]) j--;
            if (i >= j) {
                break;
            } // No more empty places

            temp = *(QDot + i);
            *(QDot + i) = *(QDot + j);
            *(QDot + j) = temp;

            // Update the arrays 
            alive[i] = true;
            alive[j] = false;
            e_local_old[i] = e_local_old[j];
            sorted++;
            j--;
        }

        i++;
        if (i >= j) break;
    }
    
    num_alive = num_alive - killed + num_copies; // Update num_alive

    return;
    
}

