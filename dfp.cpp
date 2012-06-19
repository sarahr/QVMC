/* 
 * File:   dfp.cpp
 * Author: sarahrei
 *  
 * Created on 14. march 2012, 16:42
 * 
 * Based on the DFP algorithm, the local energy is minimized with respect to 
 * the variational parameters alpha and beta. 
 * Note: Set MINIMIZE to "true" in qvmc.cpp !!!
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

/*
 * Compute the local energy based on the VMC algorithm. At the same time,
 * the gradient of the wavefunction with respect to the variational parameters 
 * is computed and stored in the vectors total_par_psi and total_par_psi2. 
 * These two vectors are then transferred to the function delE_func().
 */
double E_func(vec x, int sampling, Metropolis* VMC_brute, Metropolis_Hastings*
        VMC_imp, int N, int numprocs, int myrank, double del_min, double del_max,
        double eps, double N_therm, double N_delta, double step, int numpart,
        vec& total_par_psi, vec& total_par_psi2);

/*
 * Compute the gradient of the local energy with respect to the variational 
 * parameters. Based on the two input vectors that were computed in the 
 * function E_func().
 */
vec delE_func(vec& total_par_psi, vec& total_par_psi2);

void lnsrch(int n, vec &xold, double fold, vec &g, vec &p, vec &x,
        double *f, double stpmax, int *check, double (*func)(vec x, int sampling,
        Metropolis* VMC_brute, Metropolis_Hastings* VMC_imp, int N, int numprocs,
        int myrank, double del_min, double del_max, double eps, double N_therm,
        double N_delta, double step, int numpart, vec& total_par_psi, vec&
        total_par_psi2), int sampling, Metropolis* VMC_brute, Metropolis_Hastings*
        VMC_imp, int N, int numprocs, int myrank, double del_min, double del_max,
        double eps, double N_therm, double N_delta, double step, int numpart,
        vec& total_par_psi, vec& total_par_psi2, vec start);


void dfpmin(vec &p, int n, double gtol, int *iter, double *fret, double(*func)
        (vec x, int sampling, Metropolis* VMC_brute, Metropolis_Hastings* VMC_imp,
        int N, int numprocs, int myrank, double del_min, double del_max, double eps,
        double N_therm, double N_delta, double step, int numpart, vec& total_par_psi,
        vec& total_par_psi2), vec(*dfunc)(vec& p, vec& q), int sampling, Metropolis*
        VMC_brute, Metropolis_Hastings* VMC_imp, int N, int numprocs, int myrank,
        double del_min, double del_max, double eps, double N_therm, double N_delta,
        double step, int numpart, vec& total_par_psi, vec& total_par_psi2);


static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define ALF 1.0e-4



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
    }
    outfilename = argv[1];
    ofile.open(outfilename);


    //MPI initialization 
    int numprocs;
    int myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    int N, N_therm, N_delta; // MC cycles, Thermalization, Determine delta
    double alpha, beta; // Variational parameters
    double local_sum, local_squaresum, local_r, local_rsq;
    long idum = -time(NULL) - myrank - 1; // Seed for random number generator
    double del_max = 3.0; // Parameters to determine the optimal value
    double del_min = 0.01; // of delta
    double eps = .001; // Tolerance for delta
    vec total_par_psi(2), total_par_psi2(2);
    int numpart, dim;
    double omega, step;
    double a_start, a_delta, b_start, b_delta;
    int code; // 0 for analytical, 1 for numerical expressions
    int sampling; // 0 for brute force, 1 for importance sampling 
    int interaction; // 0: with interaction, 1: without interaction
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
    a_delta = INIreader.GetDouble("main", "a_delta");
    b_delta = INIreader.GetDouble("main", "b_delta");

    code = INIreader.GetInt("main", "code1");
    sampling = INIreader.GetInt("main", "sampling");
    interaction = INIreader.GetInt("main", "interaction");
    jastrow = INIreader.GetInt("main", "jastrow");

    // Initializations
    Wavefunction* Trial = new Wavefunction(numpart, dim, omega, jastrow);
    Hamiltonian Hamilt(numpart, dim, interaction, code, omega);
    Metropolis VMC_brute(Trial, &Hamilt, idum);
    Metropolis_Hastings VMC_imp(Trial, &Hamilt, idum);

    int nn, iter;
    double gtol, fret;

    nn = 2;
    //   reserve space in memory for vectors containing the variational
    //   parameters
    vec p(nn);

    gtol = 1.0e-5;

    p(0) = a_start;
    p(1) = b_start;
    dfpmin(p, nn, gtol, &iter, &fret, E_func, delE_func, sampling, &VMC_brute,
            &VMC_imp, N, numprocs, myrank, del_min, del_max, eps, N_therm,
            N_delta, step, numpart, total_par_psi, total_par_psi2);

    if (myrank == 0) {
        cout << "Value of energy minimum = " << fret << endl;
        cout << "Number of iterations = " << iter << endl;
        cout << "Value of alpha at minimum = " << p(0) << endl;
        cout << "Value of beta at minimum = " << p(1) << endl;

        ofile.close();
    }

    MPI_Finalize();



    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                  Functions
///////////////////////////////////////////////////////////////////////////////

double E_func(vec x, int sampling, Metropolis* VMC_brute, Metropolis_Hastings*
        VMC_imp, int N, int numprocs, int myrank, double del_min, double del_max,
        double eps, double N_therm, double N_delta, double step, int numpart,
        vec& total_par_psi, vec& total_par_psi2) {

    double local_sum = 0;
    double total_sum = 0;
    vec local_par_psi(2), local_par_psi2(2);

    total_par_psi.zeros();
    total_par_psi2.zeros();

    // Distribution of MC cycles to processors
    N /= numprocs;

    //Run the MC algorithm
    if (sampling == 0) {
        VMC_brute->delta_opt(N_delta, x(0), x(1), del_min, del_max, eps);
        VMC_brute->run_algo(N, N_therm, x(0), x(1), myrank);
        local_sum = VMC_brute->get_E();
        local_par_psi = VMC_brute->get_par_psi();
        local_par_psi2 = VMC_brute->get_par_psi2();

    } else {
        VMC_imp->set_delt(step);
        VMC_imp->run_algo(N, N_therm, x(0), x(1), myrank);
        local_sum = VMC_imp->get_E();
        local_par_psi = VMC_imp->get_par_psi();
        local_par_psi2 = VMC_imp->get_par_psi2();
    }


    // Collect data from all processors
    MPI_Allreduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int k = 0; k < 2; k++) {
        MPI_Allreduce(&local_par_psi(k), &total_par_psi(k), 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_par_psi2(k), &total_par_psi2(k), 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
    }


    N *= numprocs;
    double E = total_sum / (N * numpart);
    total_par_psi /= N * numpart;
    total_par_psi2 /= N * numpart;
    total_par_psi *= E;

    return E;

}

vec delE_func(vec& total_par_psi, vec& total_par_psi2) {

    vec ret(2);
    ret = 2 * total_par_psi2 - 2 * total_par_psi;

    return ret;
}

void dfpmin(vec &p, int n, double gtol, int *iter, double *fret, double(*func)
        (vec x, int sampling, Metropolis* VMC_brute, Metropolis_Hastings* VMC_imp,
        int N, int numprocs, int myrank, double del_min, double del_max, double eps,
        double N_therm, double N_delta, double step, int numpart, vec& total_par_psi,
        vec& total_par_psi2), vec(*dfunc)(vec& x, vec& y), int sampling, Metropolis*
        VMC_brute, Metropolis_Hastings* VMC_imp, int N, int numprocs, int myrank,
        double del_min, double del_max, double eps, double N_therm, double N_delta,
        double step, int numpart, vec& total_par_psi, vec& total_par_psi2) {

    double a_start = p(0);
    double b_start = p(1);
    vec start(2);
    start = p;

    int check, i, its, j;
    double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test, dum;
    vec dg(n), g(n), hdg(n), pnew(n), xi(n);
    mat hessian(n, n);

    fp = (*func)(p, sampling, VMC_brute, VMC_imp, N, numprocs, myrank, del_min,
            del_max, eps, N_therm, N_delta, step, numpart, total_par_psi, total_par_psi2);
    g = (*dfunc)(total_par_psi, total_par_psi2);

    g /= norm(g, 2); // Normalize the first gradient

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) hessian(i, j) = 0.0;
        hessian(i, i) = 1.0;
        xi(i) = -g(i);
        sum += p(i) * p(i);
    }

    stpmax = STPMX * FMAX(sqrt(sum), (double) n);
    for (its = 1; its <= ITMAX; its++) {
        if (myrank == 0) cout << p << endl;
        *iter = its;
        if (p(0) < 0) p(0) = a_start; // if alpha or beta negative, take again the
        if (p(1) < 0) p(1) = b_start; // starting parameters

        lnsrch(n, p, fp, g, xi, pnew, fret, stpmax, &check, func,
                sampling, VMC_brute, VMC_imp, N, numprocs, myrank, del_min,
                del_max, eps, N_therm, N_delta, step, numpart, total_par_psi,
                total_par_psi2, start);
        fp = *fret;
        for (i = 0; i < n; i++) {
            xi(i) = pnew(i) - p(i);
            p(i) = pnew(i);
        }


        test = 0.0;
        for (i = 0; i < n; i++) {
            temp = fabs(xi(i)) / FMAX(fabs(p(i)), 1.0);
            if (temp > test) test = temp;
        }
        if (test < TOLX) {
            return;
        }
        for (i = 0; i < n; i++) dg(i) = g(i);

        dum = (*func)(p, sampling, VMC_brute, VMC_imp, N, numprocs, myrank,
                del_min, del_max, eps, N_therm, N_delta, step, numpart,
                total_par_psi, total_par_psi2);

        g = (*dfunc)(total_par_psi, total_par_psi2);

        test = 0.0;
        den = FMAX(*fret, 1.0);
        for (i = 0; i < n; i++) {
            temp = fabs(g(i)) * FMAX(fabs(p(i)), 1.0) / den;
            if (temp > test) test = temp;
        }
        if (test < gtol) {
            return;
        }
        for (i = 0; i < n; i++) dg(i) = g(i) - dg(i);
        for (i = 0; i < n; i++) {
            hdg(i) = 0.0;
            for (j = 0; j < n; j++) hdg(i) += hessian(i, j) * dg(j);
        }
        fac = fae = sumdg = sumxi = 0.0;
        for (i = 0; i < n; i++) {
            fac += dg(i) * xi(i);
            fae += dg(i) * hdg(i);
            sumdg += SQR(dg(i));
            sumxi += SQR(xi(i));
        }
        if (fac * fac > EPS * sumdg * sumxi) {
            fac = 1.0 / fac;
            fad = 1.0 / fae;
            for (i = 0; i < n; i++) dg(i) = fac * xi(i) - fad * hdg(i);
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    hessian(i, j) += fac * xi(i) * xi(j)
                            - fad * hdg(i) * hdg(j) + fae * dg(i) * dg(j);
                }
            }
        }
        for (i = 0; i < n; i++) {
            xi(i) = 0.0;
            for (j = 0; j < n; j++) xi(i) -= hessian(i, j) * g(j);
        }
    }
    cout << "too many iterations in dfpmin" << endl;


}

void lnsrch(int n, vec &xold, double fold, vec &g, vec &p, vec &x,
        double *f, double stpmax, int *check, double (*func)(vec x, int sampling,
        Metropolis* VMC_brute, Metropolis_Hastings* VMC_imp, int N, int numprocs,
        int myrank, double del_min, double del_max, double eps, double N_therm,
        double N_delta, double step, int numpart, vec& total_par_psi, vec&
        total_par_psi2), int sampling, Metropolis* VMC_brute, Metropolis_Hastings*
        VMC_imp, int N, int numprocs, int myrank, double del_min, double del_max,
        double eps, double N_therm, double N_delta, double step, int numpart,
        vec& total_par_psi, vec& total_par_psi2, vec start) {

    int i;
    double a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope, sum, temp,
            test, tmplam;
    vec total_par_psi_test(2), total_par_psi2_test(2);

    *check = 0;
    for (sum = 0.0, i = 0; i < n; i++) sum += p(i) * p(i);
    sum = sqrt(sum);
    if (sum > stpmax)
        for (i = 0; i < n; i++) p(i) *= stpmax / sum;
    for (slope = 0.0, i = 0; i < n; i++)
        slope += g(i) * p(i);
    test = 0.0;
    for (i = 0; i < n; i++) {
        temp = fabs(p(i)) / FMAX(fabs(xold(i)), 1.0);
        if (temp > test) test = temp;
    }
    alamin = TOLX / test;
    alam = 1.0;
    for (;;) {
        for (i = 0; i < n; i++) x(i) = xold(i) + alam * p(i);
        for (int j = 0; j < 2; j++) {
            if (x(j) < 0) x(j) = start(j); // No negative alpha, beta !
        }

        *f = (*func)(x, sampling, VMC_brute, VMC_imp, N, numprocs, myrank, del_min,
                del_max, eps, N_therm, N_delta, step, numpart, total_par_psi_test,
                total_par_psi2_test);

        if (alam < alamin) {
            for (i = 0; i < n; i++) x(i) = xold(i);
            *check = 1;
            return;
        } else if (*f <= fold + ALF * alam * slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope / (2.0 * (*f - fold - slope));
            else {
                rhs1 = *f - fold - alam*slope;
                rhs2 = f2 - fold2 - alam2*slope;
                a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                if (a == 0.0) tmplam = -slope / (2.0 * b);
                else {
                    disc = b * b - 3.0 * a*slope;
                    if (disc < 0.0) cout << "Roundoff problem in lnsrch." << endl;
                    else tmplam = (-b + sqrt(disc)) / (3.0 * a);
                }
                if (tmplam > 0.5 * alam)
                    tmplam = 0.5 * alam;
            }
        }
        alam2 = alam;
        f2 = *f;
        fold2 = fold;
        alam = FMAX(tmplam, 0.1 * alam);
    }
}
#undef ALF
#undef TOLX
#undef ITMAX
#undef EPS
#undef STPMX



