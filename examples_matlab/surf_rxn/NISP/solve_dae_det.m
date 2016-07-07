function [X T] = solve_dae_det(paramb, Tf)
% Example that solves deterministic 3-equation model for heterogeneous surface
% reaction involving a monomer, dimer, and inert species adsorbing
% onto a surface out of gas phase. This model mimics some aspects
% of CO oxidation. 
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
%    Tf      --- final time
%    paramb  --- coefficient appearing in the system rhs   
%                Note: This is the coefficient we want to vary, and 
%                we expect Hopf bifurcation for b in [20.3 , 21.2].
%
% output: 
%
%    X --- matrix with each column the time evolution of the respective
%          system species (X(:,1) = u, X(:, 2) = v, X(:, 3) = w)
%
%    T --- the time mesh

global b;
b = paramb; 


% Initial conditions of zero coverage (based on Makeev:2002)
u0 = 0.0;
v0 = 0.0;
w0 = 0.0;
z0 = 1.0 - u0 - v0 - w0;

% initial time
tym = 0;

% time step
dt = 1e-2;

t = [0 : dt : Tf];
N = length(t);
X = zeros(N, 3);
T = zeros(N, 1);

% initialize
step = 1;
T(step) = 0.0;
X(step, 1) = u0;
X(step, 2) = v0;
X(step, 3) = w0;

while tym < Tf
    [u v w] = forwardFunction(dt, X(step, 1), X(step,2), X(step, 3));

    tym = tym + dt;

    step = step + 1;
    X(step, 1) = u;
    X(step, 2) = v;
    X(step, 3) = w;
    T(step) = tym;

end

% ----------------------------------
function [u v w] = forwardFunction(dt, u0, v0, w0)

z0 = 1.0 - u0 - v0 - w0;

% Integrate with 2nd order Runge Kutta

% Compute right hand sides
[dudt, dvdt, dwdt] = getRHS(u0, v0, w0, z0);

% Advance to mid-point
u = u0 + 0.5 * dt * dudt;
v = v0 + 0.5 * dt * dvdt;
w = w0 + 0.5 * dt * dwdt;
z = 1.0 - u - v - w;

% Compute right hand sides
[dudt dvdt dwdt] = getRHS(u, v, w, z);

% Advance to next time step
u = u0 + dt * dudt;
v = v0 + dt * dvdt;
w = w0 + dt * dwdt;
z = 1.0 - u - v - w;

% ----------------------------------
function [dudt dvdt dwdt] = getRHS(u, v, w, z)
global b 

  a = 1.6;
  c = 0.04;
  d = 1.0;
  e = 0.36;
  f = 0.016;

  dudt = a * z - c * u - 4.0 * d * u * v;
  dvdt = 2.0 * b * z^2 - 4.0 * d * u * v;
  dwdt = e * z - f * w;

