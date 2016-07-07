function logf = uq_log(pcdata, f)
% function: uq_exp.m
% purpose:  computes galerkin log of a PCE
%
% usage: expf = uq_log(pcdata, f)
%
% input:
%    
%    pcdata: struct containing basic PC parameters 
%    f  : is the input PCE
%
% output:
%    logf  : is the natural log of f
%
% Licensing:
%
%                The UQ Toolkit (UQTk) version 2.1.1
%                    Copyright (2013) Sandia Corporation
%                      http://www.sandia.gov/UQToolkit/
%
%     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
%     with Sandia Corporation, the U.S. Government retains certain rights in this software.
%
%     This file is part of The UQ Toolkit (UQTk)
%
%     UQTk is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     UQTk is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public License
%     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
%
%     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
%     Sandia National Laboratories, Livermore, CA, USA



% The idea here is to solve df/dx=f using forward Euler
% method.

% u is the initial guess of f
x = zeros(size(f));
u = x;
x(1) = f(1);
x0 = x;
u(1) = log(x(1));

% Number of integration steps
nsteps = 2^8;

% Step size
dx = (f-x) / nsteps;

for i=1:nsteps
    dxox = uq_div(pcdata,dx,x);
    u = u + dxox;
    x = x0 + i*dx;
end

logf=u;
