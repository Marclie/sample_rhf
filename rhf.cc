#include "rhf.h"

#include <memory>

using namespace std;

namespace psi::sample_rhf {

    MyRHF::MyRHF(psi::SharedWavefunction &ref_wfn, psi::Options &options) :
                                                                   // call constructor of psi4 Wavefunction class
                                                                   Wavefunction(ref_wfn, options),
                                                                   ref_wfn_(ref_wfn), // store reference wavefunction
                                                                   options_(options) // store options
    {
        // build molecule
        build_mol();

        // build one-electron integrals
        build_oei();

        // build the auxiliary basis set
        build_aux();
    }

    void MyRHF::build_mol() {
        // grab molecule from reference wavefunction
        SharedMolecule mol = ref_wfn_->molecule();

        // get number of atoms, electrons, and occupied orbitals
        charge_ = mol->molecular_charge(); // charge of molecule
        nelec_ = 0; // initialize number of electrons
        for (int i = 0; i < mol->natom(); i++) {
            nelec_ += (int) mol->Z(i); // add atomic number of atoms i to number of electrons
        }
        nelec_ -= charge_; // subtract charge of molecule from the number of electrons
        nocc_ = nelec_ / 2; // the number of occupied orbitals is half the number of electrons (closed-shell)

        // check that number of electrons is even
        if (nelec_ % 2 != 0) throw PsiException("Number of electrons is not even.", __FILE__, __LINE__);
        multiplicity_ = 1; // multiplicity must be 1 for closed-shell molecules

        // get nuclear repulsion energy
        array<double, 3> dipole_field = {0.0,0.0,0.0}; // dipole field (we don't need this, so it's just a placeholder)
        enuc_ = ref_wfn_->molecule()->nuclear_repulsion_energy(dipole_field);

    }

    void MyRHF::build_oei() {

        // get primary basis set from reference wavefunction
        auto basis = ref_wfn_->get_basisset("ORBITAL");

        // get number of basis functions
        nso_ = basis->nbf();

        // get one-electron integrals
        auto mints = ref_wfn_->mintshelper(); // get MintsHelper object (molecular integral helper)
        T_ = mints->ao_kinetic(); // kinetic energy integrals in AO basis (atomic orbital basis)
        V_ = mints->ao_potential(); // potential energy integrals in AO basis (atomic orbital basis)
        S_ = mints->ao_overlap(); // overlap integrals in AO basis (atomic orbital basis)

        // build core Hamiltonian
        H_ = T_->clone(); // initialize core Hamiltonian as kinetic energy integrals (make a copy of T_)
        H_->add(V_); // add potential energy integrals to core Hamiltonian (H_ = T_ + V_)

        build_X(); // build orthogonalizer
    }

    void MyRHF::build_X(){
        /** build orthogonalizer that rotates AO basis to an orthonormal basis (called the SO basis)
         * this is a very important step in SCF
         * the AO basis is not orthonormal, so we need to rotate it to an orthonormal basis.
         * this is done by diagonalizing the overlap matrix (S), taking the inverse square root of each eigenvalue,
         * and then transforming the new eigenvalues back to the original basis:
         *      s = U^T * S * U                 <- s is a diagonal matrix of eigenvalues
         *      X = U * s^(-1/2) * U^T
         * X is the orthogonalizer such that
         *      X^T * S * X = I (identity matrix)
         * */
        X_ = S_->clone(); // initialize orthogonalizer as overlap integrals (make a copy of S_)
        X_->power(-0.5); // diagonalize overlap integrals, inverse square root each eigenvalue, transform back
    }

    void MyRHF::build_aux() {

        // get primary basis set from reference wavefunction
        auto basis = ref_wfn_->get_basisset("ORBITAL");

        // get auxiliary basis set
        auto auxbasis = ref_wfn_->get_basisset("DF_BASIS_SCF");

        // get number of auxiliary basis functions
        nQ_ = auxbasis->nbf();

        // initialize C matrix
        C_ = make_shared<Matrix>(nso_, nso_);

        // build density-fitting tensor
        int nvir = nso_ - nocc_; // number of unoccupied orbitals (virtual orbitals)
        shared_ptr<DFTensor> DF (new DFTensor(basis, auxbasis, C_, nocc_, nvir, nocc_, nvir, options_));
        Qso_ = DF->Qso(); // get three-index tensor: (μν|Q) and (Q|λσ) integrals, where Qso is symmetric
    }

    double MyRHF::compute_energy() {

        // grab convergence parameters from Psi4 Options
        e_conv_ = options_.get_double("E_CONVERGENCE"); // energy convergence criteria
        d_conv_ = options_.get_double("D_CONVERGENCE"); // density convergence criteria
        maxiter_ = options_.get_int("MAXITER"); // maximum number of iterations

        // build guess wavefunction and energy
        energy_ = build_wfn();

        // print initial information
        outfile->Printf("\n    Guess Energy: %20.14f\n", energy_);

        // print header for table of iterations
        outfile->Printf("\n    ====> Begin RHF Iterations <====\n");
        outfile->Printf("\n    %5s %20s %20s %20s\n", "Iter", "Energy", "dE", "dD");

        // loop until convergence or maximum iterations is reached
        last_energy_ = energy_; // initialize last energy to zero
        last_D_ = D_->clone(); // initialize last density matrix to zero

        bool converged = false; // initialize convergence flag to false
        iter_ = 0; // initialize iteration counter for SCF loop
        do { // while not converged

            // build next wavefunction and compute energy (this is the SCF step)
            energy_ = build_wfn();

            // convergence criteria
            converged = check_convergence();

            // print iteration information
            outfile->Printf("    %5d %20.14f %20.14f %20.14f\n", iter_, energy_, delE_, delD_);

            // increment iteration counter and check if we have exceeded the maximum iterations
            if (iter_++ == maxiter_) throw PsiException("Maximum number of iterations exceeded.", __FILE__, __LINE__);
        } while (!converged); // end SCF loop

        // print footer for table of iterations
        outfile->Printf("\n    ====> End RHF Iterations <====\n");
        outfile->Printf("\n    Final RHF Energy: %20.14f\n", energy_);

        // set energy in Psi4 program
        psi::Process::environment.globals["CURRENT ENERGY"] = energy_;

        return energy_;
    }

    bool MyRHF::check_convergence() {

        /* change in energy */
        delE_ = energy_ - last_energy_; // compute the magnitude of the change in energy
        last_energy_ = energy_; // store last energy

        /* change in density */
        SharedMatrix D_diff = D_->clone(); // copy density matrix to temporary matrix
        D_diff->subtract(last_D_); // subtract last density matrix
        delD_ = D_diff->rms(); // compute root-mean-square deviation of change in density matrix
        last_D_ = D_->clone(); // copy density matrix to last density matrix

        // return true if change in energy and change in density are below convergence criteria
        return fabs(delE_) < e_conv_ && fabs(delD_) < d_conv_;

    }

    double MyRHF::build_wfn() {

        // build wavefunction quantities
        build_F(); // build Fock matrix
        build_C(); // build matrix of AO->MO coefficients
        build_D(); // build density matrix

        // compute energy
        /** the energy is given by
         *      E = D[μ,ν] (H[μ,ν] + F[μ,ν]) + E_nuc
         * where D is the density matrix, H is the core Hamiltonian, F is the Fock matrix,
         *       and E_nuc is the nuclear repulsion energy
         * */

        // initialize energy with nuclear repulsion energy
        double energy = enuc_; // initialize energy with nuclear repulsion energy

        // use vector_dot to compute the sum of the element-wise product of two matrices (D and H or F in this case)
        energy += D_->vector_dot(H_); // add contribution from core Hamiltonian
        energy += D_->vector_dot(F_); // add contribution from Fock matrix

        // return energy
        return energy;
    }

    void MyRHF::build_F() {
        /** build Fock matrix in the AO basis
         *      F = H + 2J - K
         * where H is the core Hamiltonian, J is the coulomb matrix, and K is the exchange matrix
         * */

        // check whether the guess has already been built (static, so it exists between function calls)
        static bool built_guess = false;
        if (built_guess) { // if this is not the first iteration, we need to add two-electron integrals
            // build two-electron integrals (J and K)
            build_JK();

            // Add J and K to form Fock matrix
            F_ = J_->clone(); F_->scale(2); // initialize Fock matrix as coulomb matrix and scale by 2
            F_->subtract(K_); // subtract exchange matrix

            // add core Hamiltonian to Fock matrix
            F_->add(H_);
        } else {
            F_ = H_->clone(); // else, initialize Fock matrix as core Hamiltonian
            built_guess = true; // set built_guess to true
        }
    }

    void MyRHF::build_C() {

        /** this is a generalized eigenvalue problem of the form
         *       F C = ε S C
         * where F is the Fock matrix, C is the matrix of AO->MO coefficients,
         * ε is the diagonal matrix of orbital energies, and S is the overlap matrix.
         *
         * the orthogonalizer, X, is used to transform the problem into a standard eigenvalue problem of the form
         *      F' C' = ε C'
         * where F' = X^T F X
         * and   C' = X^T C     <- this is the matrix of SO->MO coefficients
         * */

        // transform Fock matrix to the SO basis
        SharedMatrix Fs = F_->clone(); // initialize Fock matrix in the AO basis
        Fs->transform(X_); // transform Fock matrix to the SO basis

        // diagonalize Fock matrix in the SO basis
        SharedMatrix Cs = make_shared<Matrix>(nso_, nso_); // initialize matrix of MO coefficients (eigenvectors of Fock matrix)
        SharedVector eps = make_shared<Vector>(nso_); // initialize vector of orbital energies (eigenvalues of Fock matrix)
        Fs->diagonalize(Cs, eps); // diagonalize Fock matrix in the SO basis

        // transform MO coefficients back to the AO basis
        C_ = Cs->clone(); // initialize matrix of AO->MO coefficients as the matrix of SO->MO coefficients

        // gemm is a BLAS function that performs general matrix-matrix multiplication
        //      C = α A B + β C
        C_->gemm(false, false, 1.0, X_, Cs, 0.0);
    }

    void MyRHF::build_D() {
        /** the density matrix is given by
         *      D[μ,ν] = C[μ,i] C[ν,i]
         * where C is the matrix of AO->MO coefficients, and i is the number of occupied orbitals
         * */

        // initialize density matrix
        D_ = make_shared<Matrix>(nso_, nso_);
        double **D_p = D_->pointer(); // pointer to density matrix
        double **C_p = C_->pointer(); // pointer to matrix of AO->MO coefficients

        // loop over basis functions with index μ
        for (int mu = 0; mu < nso_; ++mu) {
            // loop over basis functions with index ν
            for (int nu = 0; nu < nso_; ++nu) {
                D_p[mu][nu] = 0.0; // initialize density matrix to zero

                // loop over occupied orbitals with index i
                for (int i = 0; i < nocc_; ++i) {
                    D_p[mu][nu] += C_p[mu][i] * C_p[nu][i]; // add contribution to density matrix element
                } // end loop over i
            } // end loop over ν
        } // end loop over μ
    }

    void MyRHF::build_JK() {

        /** We use the following equations to build the two-electron integrals:
         *      J[μ,ν] = D[λ,σ] * (λσ|Q) * (μν|Q)
         *      K[μ,ν] = D[λ,σ] * (μσ|Q) * (λν|Q)
         *  This is represented in einsum notation, which is a compact way of writing tensor contractions:
         *      C[i,j] = A[i,k] * B[k,j]                  ->  each i,j index of C is a sum over the κ index of A and B
         *      C[i,j] = sum_{k}( A[i,k] * B[k,j] )       ->  indices that do not appear in the output are summed over
         *
         *  We can improve the efficiency of this calculation by creating intermediate tensors:
         *      Jtmp[Q] = D[λ,σ] * (Q|λσ)
         *      Ktmp[Q,ν,σ] = D[λ,σ] * (Q|λν)
         *
         *  now we can write the two-electron integrals as:
         *    J[μ,ν] = Jtmp[Q] * (μν|Q)
         *    K[μ,ν] = Ktmp[Q,ν,σ] * (λσ|Q)
         * */

        /// build coulomb matrix
        build_J();

        /// build exchange matrix
        build_K();

    }

    void MyRHF::build_J() {

        // grab pointers to data
        double **D_p = D_->pointer(); // pointer to density matrix
        double **Qso_p = Qso_->pointer(); // pointer to Qso tensor

        /// J[μ,ν] = Jtmp[Q] * (μν|Q)

        J_ = make_shared<Matrix>(nso_, nso_); // initialize coulomb matrix
        SharedVector tmpJ(new Vector(nQ_)); // initialize intermediate for coulomb matrix

        double **J_p = J_->pointer(); // pointer to coulomb matrix
        double *tmpJ_p = tmpJ->pointer(); // pointer to intermediate for coulomb matrix

        // build intermediate for coulomb matrix

        for (int Q = 0; Q < nQ_; Q++) { // loop over auxiliary basis functions with index Q, from 0 to nQ_
            // initialize Jtmp[Q] to zero (this is a vector, so we use the pointer to the data)
            tmpJ_p[Q] = 0.0;

            for (int lam = 0; lam < nso_; lam++) { // loop over basis functions with index λ
                for (int sig = 0; sig < nso_; sig++) { // loop over basis functions with index σ

                    /** get compound index for Qso tensor
                     * Qso tensor is a three-index tensor, so we need to use compound indexing, where
                     * multiple indices are represented by a single index via the formula:
                     *    idx = idx1 + (idx2 * n1)  + (idx3 * n1*n2) + (idx4 * n1*n2*n3) + ...
                     * where n1, n2, n3, ... are the dimensions of each index
                     * */
                    int lamsig = (lam * nso_) + sig; // get compound index for Qso tensor
                    tmpJ_p[Q] += D_p[lam][sig] * Qso_p[Q][lamsig]; // add contribution to Jtmp[Q]

                } // end loop over λ
            } // end loop over σ

        } // end loop over Q
        // end build intermediate for coulomb matrix

        // build coulomb matrix

        for (int mu = 0; mu < nso_; mu++) { // loop over basis functions with index μ
            for (int nu = 0; nu < nso_; nu++) { // loop over basis functions with index ν
                // initialize J[μ,ν] to zero (this is a matrix, so we use the pointer to the data)
                J_p[mu][nu] = 0.0;

                for (int Q = 0; Q < nQ_; Q++) { // loop over auxiliary basis functions with index Q
                    int munu = (mu * nso_) + nu; // get compound index for J matrix
                    J_p[mu][nu] += tmpJ_p[Q] * Qso_p[Q][munu]; // add contribution to J[μ,ν]
                } // end loop over Q

            } // end loop over ν
        } // end loop over μ
        // end build coulomb matrix

    }

    void MyRHF::build_K() {

        // grab pointers to data
        double **D_p = D_->pointer(); // pointer to density matrix
        double **Qso_p = Qso_->pointer(); // pointer to Qso tensor

        /// K[μ,ν] = Ktmp[Q,ν,σ] * (λσ|Q)

        K_ = make_shared<Matrix>(nso_, nso_); // initialize exchange matrix
        SharedMatrix tmpK(new Matrix(nQ_, nso_ * nso_)); // initialize intermediate for exchange matrix
        double **K_p = K_->pointer(); // pointer to exchange matrix
        double **tmpK_p = tmpK->pointer(); // pointer to intermediate for exchange matrix

        // build intermediate for exchange matrix
        for (int Q = 0; Q < nQ_; ++Q) { // loop over auxiliary basis functions with index Q
            for (int nu = 0; nu < nso_; ++nu) { // loop over basis functions with index ν
                for (int sig = 0; sig < nso_; ++sig) { // loop over basis functions with index σ

                    int nusig = (nu * nso_) + sig; // get compound index for intermediate
                    tmpK_p[Q][nusig] = 0.0; // initialize intermediate to zero

                    // loop over basis functions with index λ
                    for (int lam = 0; lam < nso_; ++lam) {
                        int lamnu = (lam * nso_) + nu; // get compound index for Qso tensor
                        tmpK_p[Q][nusig] += D_p[lam][sig] * Qso_p[Q][lamnu]; // add contribution to intermediate
                    } // end loop over λ

                } // end loop over σ
            } // end loop over ν
        } // end loop over Q
        // end build intermediate for exchange matrix

        // build exchange matrix

        for (int mu = 0; mu < nso_; ++mu) { // loop over basis functions with index μ
            for (int nu = 0; nu < nso_; ++nu) { // loop over basis functions with index ν

                // initialize K[μ,ν] to zero (this is a matrix, so we use the pointer to the data)
                K_p[mu][nu] = 0.0;

                for (int Q = 0; Q < nQ_; ++Q) { // loop over auxiliary basis functions with index Q
                    for (int sig = 0; sig < nso_; ++sig) { // loop over basis functions with index σ
                        int musig = (mu * nso_) + sig; // get compound index for Qso tensor
                        int nusig = (nu * nso_) + sig; // get compound index for intermediate
                        K_p[mu][nu] += tmpK_p[Q][nusig] * Qso_p[Q][musig]; // add contribution to K[μ,ν]
                    } // end loop over Q
                } // end loop over σ

            } // end loop over ν
        } // end loop over μ
        // end build exchange matrix
    }

} // End namespace
