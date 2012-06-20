/* 
 * File:   main_DMC.cpp
 * Author: sarahrei
 * 
 * Created on April 23, 2012, 7:25 PM
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "DMC.h"
#include "ini.h"
#include <time.h>
#include <mpi.h>
#include "normal.hpp"

#define PAIRCOR true // set "true" for pair correlation function

using namespace std;
using namespace arma;


/////////////////////////////////////////////////////////////////////////       
//                         The main program 
/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    int numprocs;
    int myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int N_therm = 10000;
    int N = atoi(argv[1]);
    int N_equi = atoi(argv[2]);
    double alpha = atof(argv[3]);
    double beta = atof(argv[4]);
    double omega = atof(argv[5]);
    int nw_init = atoi(argv[6]);
    double E_t = atof(argv[7]);

    double del_t = 0.002;
    int numpart = atoi(argv[8]);
    int samples = atoi(argv[9]);
    char* output = argv[10];

    int dim = 2;
    int code = 0;
    int interaction = 1;
    int jastrow = 1;
    int equi_samples = 120;

    nw_init /= numprocs;

    // Initializations
    Hamiltonian Hamilt(numpart, dim, interaction, code, omega);
    DMC DiffMonCa(numpart, dim, omega, jastrow, &Hamilt,
            alpha, beta, myrank);
    

    DiffMonCa.run_algo(N, N_therm, del_t, E_t, samples,
            N_equi, equi_samples, nw_init, output, myrank, numprocs);

    
    MPI_Finalize();

    return 0;
}

