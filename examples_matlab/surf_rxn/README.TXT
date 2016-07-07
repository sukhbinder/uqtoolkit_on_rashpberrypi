This folder contains codes for intrusive and non-intrusive UQ
for the stochastic 3-equation model for heterogeneous surface
reaction involving a monomer, dimer, and inert species adsorbing
onto a surface out of gas phase. This model mimics some aspects
of CO oxidation. The parameter 'b' of the model below is considered random,
and intrusive spectral projectio is used to solve for time evolution of the
PC modes of the solution.
%
For more details on the reaction model, see
[1] R. Vigil and F. Willmore, “Oscillatory dynamics in a heterogeneous surface
    reaction: Breakdown of the mean-field approximation.,” Phys Rev E,
    vol. 54, no. 2, pp. 1225–1231, Aug. 1996.
[2] A. G. Makeev, D. Maroudas, and I. G. Kevrekidis, “‘Coarse’ stability and
    bifurcation analysis using stochastic simulators: Kinetic Monte Carlo examples,”
    J. Chem. Phys., vol. 116, no. 23, p. 10083, 2002.
The equations solved are:
    du/dt = az - cu - 4duv (coverage fraction of monomer)
    dv/dt = 2bz^2 - 4duv   (coverage fraction of dimer)
    dw/dt = ez - fw        (coverage fraction of inert species)
      z   = 1 - u - v - w  (vacant fraction)
 
The respective folders for ISP and NISP contain the codes for intrusive and
non-intrusive uncertainty propagation respectively.
In each case a driver routine is included to illustrate the usage.
