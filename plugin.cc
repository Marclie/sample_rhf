/*
 * @BEGIN LICENSE
 *
 * rhf by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <utility>

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libqt/qt.h"
#include "rhf.h"

namespace psi::sample_rhf {

    extern "C" PSI_API
    int read_options(std::string name, Options& options)
    {
        if (name == "SAMPLE_RHF"|| options.read_globals()) {
            /*- The amount of information printed to the output file -*/
            options.add_int("PRINT", 1);
        }

        return true;
    }

    extern "C" PSI_API
    SharedWavefunction sample_rhf(SharedWavefunction ref_wfn, Options& options)
    {
        outfile->Printf("\n\n ==> RHF <==\n\n");

        // data structure for RHF
        auto my_rhf = std::make_shared<MyRHF>(ref_wfn, options);

        // compute energy
        my_rhf->compute_energy();

        return my_rhf;
    }
} // End namespaces

