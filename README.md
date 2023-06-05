This sample Restricted Hartree-Fock Psi4 plugin with descriptive comments for use in the DePrince Research Group.
It is intended to be used as a teaching tool for learning how to write Psi4 plugins using the Psi4 C++ API 
and how to use object-oriented programming to abstract away the details of the SCF algorithm.
It is not intended to be used for any serious calculations and does not include any optimizations.

## Installation
1. Clone this repository into your Psi4 plugin directory.
2. Run `psi4 --plugin-compile | /bin/bash` to generate the plugin Makefile.
3. Run `make` to compile the plugin.

## Usage
4. Run `psi4 input.dat` to run the plugin and generate the output file.
5. Compare the energy of the plugin to the reference energy from the Psi4 output.
6. Learn a little bit more about Psi4 plugins!

Refer to the [group website](https://www.chem.fsu.edu/~deprince/programming_projects/scf/) for more information.
A tutorial video for this specific plugin is pending.