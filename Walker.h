/* 
 * File:   Walker.h
 * Author: sa_rei
 *
 * Created on April 22, 2012, 11:52 PM
 */

#ifndef WALKER_H
#define WALKER_H

#include "Radial.h"
#include "Hermite.h"
#include "Jastrow.h"
#include "Slater.h"
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "QForce.h" 
#include "zignor.h" 
#include "zigrandom.h"
#include "ziggurat.hpp"
#include "normal.hpp"

#include</mn/felt/u9/sarahrei/General/Libraries/usr/include/armadillo>
//#include <armadillo> 

using namespace std;
using namespace arma;

class Walker {
public:

    int numpart;
    int dim;
    Wavefunction* Trial;
    Hamiltonian* H;
    QForce* QFo;

    double alpha, beta;
    double del_t;
    double E_old, E_new;
    long idum;
    int seed;
    mat R_cur, R_tr, trial_inv;
    bool interaction;
    double wf_R;

public:

    friend class ControlWalkers;

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
    Walker(int numpart, int dim, double omega, int jastrow,
            Hamiltonian* H, double alpha, double beta);

    /**
     * Destructor
     */
    ~Walker();

    /**
     * Initialization of walker
     * @param del_t - time step
     */
    void initWalker(double del_t);




    /**
     * If the move is not accepted, reset position and Slater inverse
     * @param p - particle that has been moved
     */
    void not_accept(int p);

    /**
     * Computation of the trial position according to the Metropolis-Hastings algorithm
     * @param p - particle that has been moved compared to previous configuration
     */
    void trial_pos(int p);

    /**
     * Compute acceptance ratio according to the Metropolis-Hastings algorithm
     * @param p - particle that has been moved 
     * @return acceptance ratio 
     */
    double ratio(int p);

    /**
     * Performs all updates if a move has been accepted
     * @param p - particle that has been moved 
     */
    void accept(int p);


    /**
     * Compute the local energy
     * @return local energy
     */
    double E_local();


    /**
     * Check if a node has been crossed
     * @return "true" if node crossed, otherwise "false"
     */
    bool nodeCrossed();



};
#endif

