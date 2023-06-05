#ifndef RHF_H
#define RHF_H

#include <utility>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"

namespace psi::sample_rhf {
    class MyRHF : public Wavefunction {

    public:
        SharedWavefunction ref_wfn_;
        Options &options_;

        // constructor
        MyRHF(SharedWavefunction &ref_wfn, Options &options);
        double compute_energy() override; // SCF procedure

    private:

        // molecule
        int natom_; // number of atoms
        int nelec_; // number of electrons
        int nocc_; // number of occupied orbitals
        int charge_; // charge of molecule
        int multiplicity_; // multiplicity of molecule
        double enuc_; // nuclear repulsion energy
        void build_mol(); // build molecule

        // energy of molecule
        double energy_;

        // basis set
        int nso_; // number of basis functions
        int nQ_; // number of auxiliary basis functions

        // one-electron integrals
        SharedMatrix T_; // kinetic energy
        SharedMatrix V_; // potential energy
        SharedMatrix H_; // core Hamiltonian
        void build_oei(); // build one-electron integrals

        SharedMatrix S_; // overlap matrix
        SharedMatrix X_; // orthogonalizer (S^-1/2)
        void build_X(); // build orthogonalizer

        // two-electron integrals
        SharedMatrix Qso_; // three-index electron repulsion integrals
        void build_aux(); // build the auxiliary basis set

        SharedMatrix J_; // coulomb matrix
        SharedMatrix K_; // exchange matrix
        void build_JK(); // build two-electron integrals
        void build_J(); // build coulomb matrix
        void build_K(); // build exchange matrix

        // matrices that determine the wavefunction
        SharedMatrix F_; // Fock matrix
        SharedMatrix C_; // MO coefficients
        SharedMatrix D_; // density matrix
        build_F(); // build Fock matrix
        build_C(); // build MO coefficients
        build_D(); // build density matrix

        // build wavefunction quantities from current density matrix and return the energy
        double build_wfn();

        // SCF convergence
        double e_conv_; // energy convergence
        double d_conv_; // density convergence
        int maxiter_; // maximum number of iterations
        int iter_; // current iteration

        double last_energy_; // last energy
        SharedMatrix last_D_; // last density matrix
        double delE_, delD_; // change in energy/density
        bool check_convergence();
    };
} // End namespaces

#endif
