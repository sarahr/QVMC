#include "qvmc.h"

//////////////////////////////// ///////////////////////////////////////////////
//                          class Hamiltonian

/** @brief Class that unifies kinetic, potential and interaction energy of the
 * Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 * @param dim -  dimension
 * @param Kin - pointer to an object of class Kinetic
 * @param Pot - pointer to an object of class Potential
 * @param Int - pointer to an object of class Interaction
 * @param intaction - 0: without interaction, 1: with interaction
 */
Hamiltonian::Hamiltonian(int numpart, int dim, Kinetic* Kin, Potential* Pot,
        Interaction* Int, int intaction) {

    this->numpart = numpart;
    this->dim = dim;
    this->Kin = Kin;
    this->Pot = Pot;
    this->Int = Int;

    if (intaction == 1) interaction = true;
    else interaction = false;

}

/**
 * Compute the unperturbed part of the Hamiltonian
 * @param Psi - pointer to wavefunction
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return  H_0
 */
double Hamiltonian::H_0(Wavefunction* Psi, double alpha, double beta) {
    return Kin->Laplace(Psi, alpha, beta) + Pot->getPotential(Psi);

}

/**
 * Compute the interaction part of the Hamiltonian
 * @param Psi - pointer to wavefunction
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return  H_1
 */
double Hamiltonian::H_1(Wavefunction* Psi) {
    return Int->getInteraction(Psi);
}


////////////////////////////////////////////////////////////////////////////////
//                              class HarmOs

/** @brief Subclass of Potential: Harmonic oscillator potential
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param omega - osc. frequency omega
 * @param numpart - number of particles
 */
HarmOs::HarmOs(double omega, int numpart) {
    this->omega = omega;
    this->numpart = numpart;
}

/**
 * Computes the potential energy of a given wavefunction, this one is given by
 * the harmonic oscillator potential
 * @param Psi - pointer to wavefunction
 * @return potential energy
 */
double HarmOs::getPotential(Wavefunction* Psi) {

    double pot = 0.0;

    for (int i = 0; i < numpart; i++) {
        pot += Psi->Pos->r(i);
    }

    pot *= 0.5 * omega*omega;

    return pot;

}


////////////////////////////////////////////////////////////////////////////////
//                           class Interaction

/** @brief Class for the interaction part of the Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param numpart - number of particles
 */
Interaction::Interaction(int numpart) {
    this->numpart = numpart;
}

/**
 * Computes the interaction energy of a given wavefunction, this one is given
 * by the dimensionless Coulomb repulsion
 * @param Psi - pointer to wavefunction
 * @return interaction energy
 */
double Interaction::getInteraction(Wavefunction* Psi) {

    int i, j;
    double intact = 0.0;

    for (j = 1; j < numpart; j++) {
        for (i = 0; i < j; i++) {

            intact += 1.0 / Psi->Pos->r_int(i, j);
        }
    }

    return intact;

}

////////////////////////////////////////////////////////////////////////////////
//                          class Kinetic

/** @brief Class for the kinetic part of the Hamiltonian
    @author sarahrei
    @date 11 April 2012
 */
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor
 * @param code - 0: analytical, 1: numerical expression for Laplacian
 * @param dim - dimension
 * @param numpart - number of particles
 * @param omega - osc. frequency omega
 * @param intaction - 0: without interaction, 1: with interaction
 */
Kinetic::Kinetic(int code, int dim, int numpart, double omega, int intaction) {

    this->dim = dim;
    this->numpart = numpart;
    this->omega = omega;
    engloc = eloc_table[code];
    h = 0.001;
    h2 = 1000000;

    if (intaction == 1) interaction = true;
    else interaction = false;
}



/*
 * Defines an array of two functions:
 * 1st component: Analytical calculation of the Laplacian
 * 2nd component: Numerical calculation of the Laplacian
 */
const Kinetic::fptr_eloc Kinetic::eloc_table[] = {
    &Kinetic::Analytical, &Kinetic::Numerical
};

/**
 * Analytical expression for the kinetic part of the local energy.
 * Usage speeds up the code. 
 * @param Psi - Pointer to wavefunction
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return - kinetic part of the local energy
 */
double Kinetic::Analytical(Wavefunction* Psi, double alpha, double beta) {

    double e_kin = 0.0;
    double lap_jas = 0;
    double lap_slater = 0;

    for (int p = 0; p < numpart; p++) {

        lap_jas += Psi->JastrowPsi->laplace(Psi->Pos, p, beta);
        lap_slater += Psi->SlaterPsi->laplace(Psi->Pos, p, alpha);
        e_kin += Psi->SlaterPsi->laplace(Psi->Pos, p, alpha);

        // Terms that arise from the interaction
        if (interaction) e_kin += Psi->JastrowPsi->laplace(Psi->Pos, p, beta);
        if (interaction) e_kin += 2 * dot(Psi->SlaterPsi->gradient(
                Psi->Pos->current, p, alpha, Psi->SlaterPsi->slat_inv),
                Psi->JastrowPsi->gradient(Psi->Pos, p, beta));

    }
    e_kin *= -0.5;

    return e_kin;
};

/**
 * Numerical computation of the kinetic energy.
 * The second derivative is discretized with truncation error O(h^2).
 * @param Psi - Pointer to wavefunction
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return - kinetic part of the local energy
 */
double Kinetic::Numerical(Wavefunction* Psi, double alpha, double beta) {

    double kinetic = 0.0;

    mat R = Psi->Pos->current;
    mat R_plus = zeros(numpart, dim);
    mat R_minus = zeros(numpart, dim);

    double wf_minus, wf_plus, wf_old;
    R_plus = R;
    R_minus = R;

    wf_old = Psi->value(R, alpha, beta);

    for (int j = 0; j < dim; j++) {
        for (int i = 0; i < numpart; i++) {

            R_plus(i, j) = R(i, j) + h;
            R_minus(i, j) = R(i, j) - h;

            wf_minus = Psi->value(R_minus, alpha, beta);
            wf_plus = Psi->value(R_plus, alpha, beta);

            kinetic -= (wf_minus + wf_plus - 2 * wf_old);

            R_plus(i, j) = R(i, j);
            R_minus(i, j) = R(i, j);
        }
    }

    // Finalize the second derivative and divide by wave function 
    kinetic = 0.5 * h2 * kinetic / wf_old;

    wf_old = Psi->value(R, alpha, beta); // to reset internal values in object

    return kinetic;

};

/**
 * Returns the kinetic part of the local energy (predominated by the Laplacian)
 * The pointer "engloc", which has been set in the constructor, picks the right
 * computation scheme (analytically or numerically)
 * @param Psi - Pointer to wavefunction
 * @param alpha - first variational parameter
 * @param beta - second variational parameter
 * @return - kinetic part of the local energy
 */
double Kinetic::Laplace(Wavefunction* Psi, double alpha, double beta) {
    return (this->*engloc)(Psi, alpha, beta);
}
