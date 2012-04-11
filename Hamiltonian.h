/* 
 * File:   Hamiltonian.h
 * Author: sa_rei
 *
 * Created on April 9, 2012, 7:45 PM
 */

#ifndef HAMILTONIAN_H
#define	HAMILTONIAN_H

#include <armadillo> 

using namespace std;
using namespace arma;


/** @brief Class for the kinetic part of the Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
class Kinetic{  
private:
    int dim;
    int numpart;
    double h; 
    double h2; 
    double omega;
    bool interaction;
    
     // Definitions for the E_local array
    typedef double (Kinetic::*fptr_eloc)(Wavefunction*,double, double);
    static const fptr_eloc eloc_table[2];   
    
public:
    
    /**
     * Constructor
     * @param code - 0: analytical, 1: numerical expression for Laplacian
     * @param dim - dimension
     * @param numpart - number of particles
     * @param omega - osc. frequency omega
     * @param intaction - 0: without interaction, 1: with interaction
     */
    Kinetic(int code, int dim, int numpart, double omega, int intaction);
    
    fptr_eloc engloc;
    
    /**
     * Analytical expression for the kinetic part of the local energy.
     * Usage speeds up the code. 
     * @param Psi - Pointer to wavefunction
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return - kinetic part of the local energy
     */
    double Analytical(Wavefunction* Psi, double alpha, double beta);
    
    /**
     * Numerical computation of the kinetic energy.
     * The second derivative is discretized with truncation error O(h^2).
     * @param Psi - Pointer to wavefunction
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return - kinetic part of the local energy
     */
    double Numerical(Wavefunction* Psi, double alpha, double beta);
    
    /**
     * Returns the kinetic part of the local energy (predominated by the Laplacian)
     * The pointer "engloc", which has been set in the constructor, picks the right
     * computation scheme (analytically or numerically)
     * @param Psi - Pointer to wavefunction
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return - kinetic part of the local energy
     */
    double Laplace(Wavefunction* Psi, double alpha, double beta);
   
};


/** @brief Class for the interaction part of the Hamiltonian
    @author sarahrei
    @date 11  April 2012
 */
class Interaction{
private:
    int numpart;
    
public:
    
    /**
     * Constructor
     * @param numpart - number of particles
     */
    Interaction(int numpart);
    
    /**
     * Computes the interaction energy for a given wavefunction, this one is given
     * by the dimensionless Coulomb repulsion
     * @param Psi - pointer to wavefunction
     * @return interaction energy
     */
    double getInteraction(Wavefunction* Psi);
};



/** @brief Class for the potential part of the Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
class Potential{
protected:
    int numpart;
    
public:
    Potential(){};
    
    /**
     * Virtual function to compute potential
     * @param Psi - pointer to wavefunction
     * @return potential part of local energy
     */
    virtual double getPotential(Wavefunction* Psi)=0;
    
};


/** @brief Subclass of Potential: Harmonic oscillator potential
    @author sarahrei
    @date 11 April 2012
 */
class HarmOs: public Potential{
private:        
    double omega;
    
public:
    
    /**
     * Constructor
     * @param omega - osc. frequency omega
     * @param numpart - number of particles
     */
    HarmOs(double omega, int numpart); 
    
    /**
     * Computes the potential energy for a given wavefunction, this one is given by
     * the harmonic oscillator potential
     * @param Psi - pointer to wavefunction
     * @return potential energy
     */
    double getPotential(Wavefunction* Psi);
    
};


/** @brief Class that unifies kinetic, potential and interaction energy of the
 * Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
class Hamiltonian{
private:
    int numpart;
    int dim; 
    Potential* Pot;
    Kinetic* Kin;
    Interaction* Int;
    bool interaction;
    
public:
    
    /**
     * Constructor
     * @param numpart - number of particles
     * @param dim -  dimension
     * @param Kin - pointer to an object of class Kinetic
     * @param Pot - pointer to an object of class Potential
     * @param Int - pointer to an object of class Interaction
     * @param intaction - 0: without interaction, 1: with interaction
     */
    Hamiltonian(int numpart, int dim, Kinetic* Kin, Potential* Pot, Interaction*
            Int, int intaction);
    
    /**
     * Compute the unperturbed part of the Hamiltonian
     * @param Psi - pointer to wavefunction
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return  unperturbed part H_0 of the Hamiltonian
     */
    double H_0(Wavefunction* Psi, double alpha, double beta); 
    
    /**
     * Compute the interaction part of the Hamiltonian
     * @param Psi - pointer to wavefunction
     * @param alpha - first variational parameter
     * @param beta - second variational parameter
     * @return  interaction part H_1 of the Hamiltonian
     */
    double H_1(Wavefunction* Psi);                              
};


#endif	/* HAMILTONIAN_H */

