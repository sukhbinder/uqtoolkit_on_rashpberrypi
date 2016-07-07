function [U V W T pcdata] = solve_dae_ISP(Tf)
% Example that solves stochastic 3-equation model for heterogeneous surface
% reaction involving a monomer, dimer, and inert species adsorbing
% onto a surface out of gas phase. This model mimics some aspects
% of CO oxidation. The parameter 'b' of the model below is considered random, 
% and intrusive spectral projectio is used to solve for time evolution of the
% PC modes of the solution.
%
% For more details on the reaction model, see
% [1] R. Vigil and F. Willmore, “Oscillatory dynamics in a heterogeneous surface
%     reaction: Breakdown of the mean-field approximation.,” Phys Rev E,
%     vol. 54, no. 2, pp. 1225–1231, Aug. 1996.
% [2] A. G. Makeev, D. Maroudas, and I. G. Kevrekidis, “‘Coarse’ stability and
%     bifurcation analysis using stochastic simulators: Kinetic Monte Carlo examples,”
%     J. Chem. Phys., vol. 116, no. 23, p. 10083, 2002.
% The equations solved are:
%     du/dt = az - cu - 4duv (coverage fraction of monomer)
%     dv/dt = 2bz^2 - 4duv   (coverage fraction of dimer)
%     dw/dt = ez - fw        (coverage fraction of inert species)
%       z   = 1 - u - v - w  (vacant fraction)
%
% input:
%
%    Tf --- final time
%
% output:
%
%    U --- Nt x nPCTerms matrix, with each column the time evolution of
%          the PC modes of u
%    V, W  Same as above, except for v and w
%
%    T --- the time mesh


global pcdata

nord = 5;
ndim = 1;
pc_type = 'LEGENDRE';

pcdata = uq_pcset(nord, ndim, pc_type);

% Initial conditions of zero coverage (based on Makeev:2002)

u0 = zeros(1,pcdata.nPCTerms);
v0 = zeros(1,pcdata.nPCTerms);
w0 = zeros(1,pcdata.nPCTerms);
z0 = zeros(1,pcdata.nPCTerms);
z0(1) = 1.0;

% initial time
tym = 0;

% time step
dt = 1e-2;

t = [0 : dt : Tf];
N = length(t);
T = zeros(N, 1);
U = zeros(N, pcdata.nPCTerms);
V = zeros(N, pcdata.nPCTerms);
W = zeros(N, pcdata.nPCTerms);

% initialize
step = 1;
T(step) = 0.0;
U(step, :) = u0;
V(step, :) = v0;
W(step, :) = w0;

hh = waitbar(0, 'Integrating, please wait ...');
while tym < Tf
   
    % update waitbar
    waitbar(step / N);
    
    [u v w] = forwardFunction(dt, U(step, :), V(step, :), W(step, :));

    tym = tym + dt;
    step = step + 1;

    U(step, :) = u;
    V(step, :) = v;
    W(step, :) = w;
    T(step) = tym;
    
end
close(hh);

% ----------------------------------
function [u v w] = forwardFunction(dt, u0, v0, w0)
global pcdata

e1 = zeros(1, pcdata.nPCTerms);
e1(1) = 1.0;

z0 = e1 - u0 - v0 - w0;

% Integrate with 2nd order Runge Kutta

% Compute right hand sides
[dudt, dvdt, dwdt] = getRHS(u0, v0, w0, z0);

% Advance to mid-point
u = u0 + 0.5 * dt * dudt;
v = v0 + 0.5 * dt * dvdt;
w = w0 + 0.5 * dt * dwdt;
z = e1 - u - v - w;

% Compute right hand sides
[dudt dvdt dwdt] = getRHS(u, v, w, z);

% Advance to next time step
u = u0 + dt * dudt;
v = v0 + dt * dvdt;
w = w0 + dt * dwdt;
z = e1 - u - v - w;

% ----------------------------------
function [dudt dvdt dwdt] = getRHS(u, v, w, z)
global pcdata

  a = 1.6;
  c = 0.04;
  d = 1.0;
  e = 0.36;
  f = 0.016;

  % the random parameter
  b = zeros(1, pcdata.nPCTerms);
  b(1) = 20.75;
  b(2) = b(1) * 0.005;

  % Build du/dt = a*z - c*u - 4.0*d*u*v
  dudt = a*z - c*u - 4.0 * d * uq_product(pcdata, u, v);

  % Build dv/dt = 2.0*b*z*z - 4.0*d*u*v
  tmp = uq_product(pcdata, z, z);
  dvdt = 2.0 * uq_product(pcdata, b, tmp) - 4.0 * d * uq_product(pcdata, u, v);

  % Build dw/dt = e*z - f*w
  dwdt = e * z - f * w;

