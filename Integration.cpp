/*
 * Project 5: One-body density
 * Name of output file is read from screen
 */

#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <armadillo>
#include "lib.h"
#include "qvmc.h"


using namespace std;
using namespace arma;

/**
 * The function One_body() performs a "dim1 x dim2"-dimensional brute force 
 * MC integration with same integration limits for all dimensions (PDF is       
 * a uniform distribution). 
 * Note: For efficiency reasons only the sum is computed,  the integral has to 
 * be finalized outside the function (but I don't want to multiply/divide 
 * each time...)
 * 
 * Variables:
 * "dim1" : dimension of single-particle space, i.e. 3 for x,y,z
 * "dim2" : number of particles to be integrated over
 * "a" and "b" : specify the integration interval [a,b]
 * "N_MC":  Number of MC cycles
 *  mat R: Cartesian coordinates of the particles
 * @param Psi - pointer to wavefunction
 * @param dim - dimension
 * @param a - lower integration limit
 * @param b - upper integration limit (interval [a,b])
 * @param N_MC - number of MC cycles
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @param idum - seed for random number generator
 * @param numpart - number of particles
 * @param coordinates - coordinates of the fixed particle, supply matrix of 
 * dimension 0 for normalization integral
 * @return integral
 */
double One_body(Wavefunction* Psi, int dim, double a, double b,
        double N_MC, double alpha, double beta, long& idum, int numpart,
        mat coordinates);

/////////////////////////////////////////////////////////////////////////
//                         The main program 
////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    char *outfilename;
    ofstream ofile;

    // Read in output file
    if (argc <= 1) {
        cout << "Bad usage: Read in output file in same line" << endl;
        exit(1);
    }

    outfilename = argv[1];
    ofile.open(outfilename);

    //MPI initialization
    int numprocs;
    int myrank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // variables for integration
    int intpoints = 30; // number integration points in each dimension
    int N_MC = 1000000; // total number of MC cycles
    int N_local = N_MC / numprocs; // MC cycles per process
    double lim = 3; // integration interval for each variable [-lim,lim] 
    int dim = 2; // dimension

    double integral; // Monte Carlo integral
    double total_int = 0;
    double norm; // Normalization factor
    double volume; // Integration volume
    double x, y, r;
    long idum = -1 - myrank; // Seed for random number generator

    double step_length = 2.0 * lim / intpoints; // Step length
    vec coordinates(dim); // coordinates of the fixed particle

    // Initialize the parameters needed for the wavefunction
    double alpha = 0.987;
    double beta = 0.4;
    int numpart = 2;
    double omega = 1.0;
    int jastrow = 1; // 0: without, 1: with Jastrow factor
    //Christoffer Hirth is supercool!!!

    // Creation of the necessary objects
    Slater Slat(numpart, omega, dim);
    ExpFactor ExpFa(numpart, omega);
    Jastrow Jast(numpart, dim);
    Radial Pos(numpart, dim);
    Radial Pos_tr(numpart, dim);
    Wavefunction Psi(numpart, dim, &Slat, &ExpFa, &Jast, &Pos, &Pos_tr, jastrow);

    // Compute normalization integral
    mat null_m(0, 0);
    volume = pow(2 * lim, dim * numpart); // Integration volume for normalization
    norm = One_body(&Psi, dim, -lim, lim, N_MC, alpha,
            beta, idum, numpart, null_m);
    norm *= volume / N_MC;

    volume = pow(2 * lim, dim * (numpart - 1)); // Integration volume for the one-body density

    for (int i = 0; i < intpoints; i++) { // x-dimension

        x = -lim + i*step_length;
        coordinates(0) = x;

        for (int j = 0; j < intpoints; j++) { // y-dimension

            y = -lim + j *step_length;
            coordinates(1) = y;

            // Compute the integral
            total_int = 0.0;
            integral = One_body(&Psi, dim, -lim, lim, N_local, alpha, beta,
                    idum, numpart, coordinates);

            r = sqrt(x * x + y * y);

            MPI_Reduce(&integral, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            N_MC = N_local*numprocs;

            if (myrank == 0) {
                total_int *= volume / N_MC / norm;
                ofile << x << "\t " << y << "\t" << total_int << "\t" << r << endl;
            }

        }
    }


    MPI_Finalize();
    ofile.close();


    return 0;
}

/*
 * The function MC_int() performs a "dim1 x dim2"-dimensional brute force 
 * MC integration with same integration limits for all dimensions (PDF is       
 * a uniform distribution). Additionally, the variance is computed. 
 * Note: For efficiency reasons only the sum and squared sum are computed. 
 * The integral and variance have to be finalized outside the function (but 
 * I don't want to multiply/divide each time...)
 * 
 * Variables:
 * "dim1" : dimension of single-particle space, i.e. 3 for x,y,z
 * "dim2" : number of particles to be integrated over
 * "a" and "b" : specify the integration interval [a,b]
 * "N_MC":  Number of MC cycles
 *  mat R: Cartesian coordinates of the particles
 */
double One_body(Wavefunction* Psi, int dim, double a, double b,
        double N_MC, double alpha, double beta, long& idum, int numpart, mat coordinates) {


    mat R(numpart, dim);
    double contr;
    int num = coordinates.n_cols;

    //double sum_var = 0.0;
    double int_MC = 0.0;

    for (int i = 0; i < N_MC; i++) {

        for (int j = num; j < numpart; j++) {
            for (int k = 0; k < dim; k++) {
                R(j, k) = a + (b - a) * ran3(&idum); // transform variables to sample space
            }
        }

        if (num != 0) {
            for (int k = 0; k < dim; k++) {
                R(0, k) = coordinates(k);
            }
        }


        contr = pow(Psi->value(R, alpha, beta), 2);
        int_MC += contr;
        //sum_var += contr*contr;
    }

    return int_MC;

}

