InterfacingFlowstar
===================

This project is a Windows port for Flow* (http://systems.cs.colorado.edu/research/cyberphysical/taylormodels) with a DLL interface for Taylor models.

Only a subset of the functions are exported, namely:
  1. intervals, 
  2. flowpipes, 
  3. continuous systems.
This allows one to use the continuous-time verified integrator (based on Taylor models) of Flow*.

In order to compile Flow*, the following libraries are needed:
  - gmp
  - mpfr
  - gsl
  - glpk
 
The binary files were produced using 32-bit MinGW and G++ 4.8.1 with the following libs:
  - gmp 5.1.3
  - mpfr 3.1.2
  - gsl 1.9
  - glpk 4.52
