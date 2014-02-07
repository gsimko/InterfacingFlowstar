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
 
