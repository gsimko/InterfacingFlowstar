InterfacingFlowstar
===================

This project is a modified Windows port of Flow* 1.2 (http://systems.cs.colorado.edu/research/cyberphysical/taylormodels) with a DLL interface for Taylor models.

Modifications
-------------

Modifications in include.h:
  1. Added an exception class for signaling cases when no valid remainder intervals can be found.
  2. Fixed a MingW issue with regards to mkdir.

Modifications in ContinuousSystem.cpp:
  1. Added iterative interval estimation widening (and miniStep narrowing) in 
	* adaptive step size, fixed-order ContinuousSystem::reach_low_degree() method.
	* fixed step size, fixed-order ContinuousSystem::reach_non_polynomial_taylor() method.
	* adaptive stepsize , fixed-order ContinuousSystem::reach_non_polynomial_taylor() method.
  2. Added exception throwing instead of fprintf(stdout) if no valid remainder intervals can be found. (ContinuousSystem::reach_... methods)

Exported functions
------------------

The following classes are exported:
  1. intervals, 
  2. flowpipes, 
  3. continuous systems,
  4. monomials,
  5. polynomials,
  6. Taylor models,
  7. ODEs,
  8. A function for setting the cutoff threshold.

Dependencies
------------

In order to compile Flow*, the following libraries are needed:
  * gmp
  * mpfr
  * gsl
  * glpk
Pre-compiled Windows binaries can be downloaded from the releases.
