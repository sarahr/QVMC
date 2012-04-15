/* 
 * File:   main.cpp
 * Author: sarahrei
 *  
 * Created on 19. januar 2012, 12:13
 */
#include <mpi.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "qvmc.h"
#include "ini.h"
#include <time.h>


using namespace std;
using namespace arma;

#define DISTANCE false // set "true" if mean distance between two particles
// is to be computed
#define E_POT_KIN false // set "true" for analyzing E_pot & E_kin separately


/////////////////////////////////////////////////////////////////////////       
//                         The main program 
/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    char *outfilename;
    ofstream ofile;

    // Read in output file
    if (argc <= 1) {
        cout << "Bad usage: Read in output file in same line" << endl;
        exit(1);
    } else {
        outfilename = argv[1];
        ofile.open(outfilename, ios_base::app);


        //MPI initialization 
        int numprocs;
        int myrank;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


        int N, N_therm, N_delta; // MC cycles, Thermalization, Determine delta
        double alpha, beta; // Variational parameters
        double local_sum, local_squaresum;
        double local_epot, local_sq_epot, local_ekin, local_sq_ekin;
        double local_r, local_rsq;
        double total_sum, total_sq;
        double total_epot, total_sq_epot, total_ekin, total_sq_ekin;
        double total_r, total_rsq;
        int numpart, dim;
        double omega, step;
        double a_start, a_steps, a_delta, b_start, b_steps, b_delta;
        double elapsed_time;

        long idum = -time(NULL) - myrank - 1; // Seed for random number generator
        double del_max = 3.0; // Parameters to determine the optimal value
        double del_min = 0.01; // of delta
        double eps = .001; // Tolerance for delta

        int code; // 0 for analytical, 1 for numerical expressions
        int sampling; // 0 for brute force, 1 for importance sampling 
        int interaction; // 0 for "off", 1 for "on"
        int jastrow; // 0 for "off", 1 for "on"

        // Get numerical values from Ini file
        ini INIreader("QD.ini");

        N = INIreader.GetInt("main", "N");
        N_therm = INIreader.GetInt("main", "N_therm");
        N_delta = INIreader.GetInt("main", "N_delta");
        numpart = INIreader.GetInt("main", "numpart");
        dim = INIreader.GetInt("main", "dim");
        omega = INIreader.GetDouble("main", "omega");
        step = INIreader.GetDouble("main", "step");

        a_start = INIreader.GetDouble("main", "a_start");
        b_start = INIreader.GetDouble("main", "b_start");
        a_steps = INIreader.GetDouble("main", "a_steps");
        b_steps = INIreader.GetDouble("main", "b_steps");
        a_delta = INIreader.GetDouble("main", "a_delta");
        b_delta = INIreader.GetDouble("main", "b_delta");

        code = INIreader.GetInt("main", "code1");
        sampling = INIreader.GetInt("main", "sampling");
        interaction = INIreader.GetInt("main", "interaction");
        jastrow = INIreader.GetInt("main", "jastrow");

        elapsed_time = -MPI_Wtime();

        // Initializations
        Wavefunction* Trial = new Wavefunction(numpart, dim, omega, jastrow);
        Hamiltonian Hamilt(numpart, dim, interaction, code, omega);
        Metropolis VMC_brute(Trial, &Hamilt, idum);
        Metropolis_Hastings VMC_imp(Trial, &Hamilt, idum);

        for (int i = 0; i < a_steps; i++) { // Loop over alpha
            for (int j = 0; j < b_steps; j++) { // Loop over beta

                alpha = a_start + i*a_delta;
                beta = b_start + j*b_delta;

                // Distribution of MC cycles to processors
                N /= numprocs;

                //Run the MC algorithm
                if (sampling == 0) {
                    VMC_brute.delta_opt(N_delta, alpha, beta, del_min, del_max, eps);
                    VMC_brute.run_algo(N, N_therm, alpha, beta);
                    local_sum = VMC_brute.get_E();
                    local_squaresum = VMC_brute.get_Esq();
#if DISTANCE
                    local_r = VMC_brute.get_rdist();
                    local_rsq = VMC_brute.get_rdistsq();
#endif

#if E_POT_KIN
                    local_epot = VMC_brute.get_Epot();
                    local_sq_epot = VMC_brute.get_Epotsq();
                    local_ekin = VMC_brute.get_Ekin();
                    local_sq_ekin = VMC_brute.get_Ekinsq();
#endif

                } else {
                    VMC_imp.set_delt(step);
                    VMC_imp.run_algo(N, N_therm, alpha, beta);
                    local_sum = VMC_imp.get_E();
                    local_squaresum = VMC_imp.get_Esq();
#if DISTANCE
                    local_r = VMC_imp.get_rdist();
                    local_rsq = VMC_imp.get_rdistsq();
#endif
#if E_POT_KIN
                    local_epot = VMC_imp.get_Epot();
                    local_sq_epot = VMC_imp.get_Epotsq();
                    local_ekin = VMC_imp.get_Ekin();
                    local_sq_ekin = VMC_imp.get_Ekinsq();
#endif
                }

                // Collect data from all processors
                MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_squaresum, &total_sq, 1, MPI_DOUBLE, MPI_SUM, 0,
                        MPI_COMM_WORLD);
#if DISTANCE               
                MPI_Reduce(&local_r, &total_r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_rsq, &total_rsq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif  

#if E_POT_KIN
                MPI_Reduce(&local_epot, &total_epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_sq_epot, &total_sq_epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_ekin, &total_ekin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_sq_ekin, &total_sq_ekin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

                // Compute final values 
                N *= numprocs;
                double E = total_sum / (N * numpart);
                double E_sq = total_sq / (N * numpart);
                double variance = E_sq - E*E;

#if DISTANCE
                total_r /= (N * numpart);
                total_rsq /= (N * numpart);
                total_rsq -= total_r*total_r;
#endif

#if E_POT_KIN
                total_epot /= (N * numpart);
                total_sq_epot /= (N * numpart);
                total_sq_epot -= total_epot*total_epot;
                total_ekin /= (N * numpart);
                total_sq_ekin /= (N * numpart);
                total_sq_ekin -= total_ekin*total_ekin;
#endif

                elapsed_time += MPI_Wtime();

                // Print output data to file
                if (myrank == 0) {

#if DISTANCE                
                    cout << alpha << " " << beta << " " << setprecision(8) << E
                            << " " << variance << " " << total_r << " " << total_rsq << endl;
                    ofile << alpha << " " << beta << " " << setprecision(8) << E
                            << " " << variance << " " << total_r << " " << total_rsq << endl;
#elif E_POT_KIN
                    cout << alpha << " " << beta << " " << E << " " << variance
                            << " " << total_epot << " " << total_sq_epot 
                            << " " << total_ekin << " " << total_sq_ekin << endl;
                    ofile << alpha << " " << beta << " " << E << " " << variance
                            << " " << total_epot << " " << total_sq_epot 
                            << " " << total_ekin << " " << total_sq_ekin << endl;
#else

                    cout << alpha << " " << beta << " " << setprecision(8) << E
                            << " " << variance << endl;
                    cout << "Elapsed time: " << elapsed_time << endl;
                    ofile << alpha << " " << beta << " " << setprecision(8) << E
                            << " " << variance << endl;
#endif
                }

            } //End loop beta
        }//End loop alpha

        MPI_Finalize();
        ofile.close();

    }

    return 0;
}


