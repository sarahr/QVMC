/*
 * Performance of the QVMC algorithm. All the contributions to the local energy
 * are written to file for further analysis by the program "blocking_analyze".
 */

#include <mpi.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "qvmc.h"
#include "ini.h"


using namespace std;
using namespace arma;


/////////////////////////////////////////////////////////////////////////       
//                         The main program 
/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    //MPI initialization 
    int numprocs;
    int myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // String handling
    ofstream ofile;
    ostringstream ost;
    ost << "blocks_rank" << myrank << ".dat";
    ofile.open(ost.str().c_str(), ios::out | ios::binary);


    int N, N_therm, N_delta; // MC cycles, Thermalization, Determine delta
    double alpha, beta; // Variational parameters
    double local_sum, local_squaresum;
    long idum = -myrank - 1; // Seed for random number generator
    int seed = myrank + 12;
    double del_max = 3.0; // Parameters to determine the optimal value
    double del_min = 0.01; // of delta
    double eps = .001; // Tolerance for delta
    double total_sum, total_sq;
    int numpart, dim;
    double omega, step;
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

    alpha = INIreader.GetDouble("main", "a_start");
    beta = INIreader.GetDouble("main", "b_start");

    code = INIreader.GetInt("main", "code1");
    sampling = INIreader.GetInt("main", "sampling");
    interaction = INIreader.GetInt("main", "interaction");
    jastrow = INIreader.GetInt("main", "jastrow");

    // Initializations
    Slater Slat(numpart, omega, dim);
    ExpFactor ExpFa(numpart, omega);
    Jastrow Jast(numpart, dim);
    Radial Pos(numpart, dim);
    Radial Pos_tr(numpart, dim);
    Kinetic Kin(code, dim, numpart, omega, interaction);
    HarmOs HaOs(omega, numpart);
    Interaction IntAct(numpart);
    Wavefunction Trial(numpart, dim, &Slat, &ExpFa, &Jast, &Pos, &Pos_tr, jastrow);
    Hamiltonian Hamilt1(numpart, dim, &Kin, &HaOs, &IntAct, interaction);
    QForce QF(&Trial);
    Metropolis VMC_brute(&Trial, &Hamilt1, idum);
    Metropolis_Hastings VMC_imp(&Trial, &Hamilt1, seed, &QF, idum);


    // Distribution of MC cycles to processors
    N /= numprocs;
    vec all_energies(N * numpart);

    //Run the MC algorithm
    if (sampling == 0) {
        VMC_brute.delta_opt(N_delta, alpha, beta, del_min, del_max, eps);
        all_energies = VMC_brute.run_algo_blocking(N, N_therm, alpha, beta);
        local_sum = VMC_brute.get_E();
        local_squaresum = VMC_brute.get_Esq();

    } else {
        VMC_imp.set_delt(step);
        all_energies = VMC_imp.run_algo_blocking(N, N_therm, alpha, beta);
        local_sum = VMC_imp.get_E();
        local_squaresum = VMC_imp.get_Esq();
    }


    // Print output data to file
    ofile << all_energies << endl;

    MPI_Finalize();
    ofile.close();

    return 0;
}

