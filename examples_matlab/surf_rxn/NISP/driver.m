function [T U V W pcdata] = driver 
% define the path to MATLAB UQ ToolBox 
% modify as needed.
addpath('../../../src_matlab/'); 

% setup PC data
nord = 5;
ndim = 1;
pc_type = 'LEGENDRE';
pcdata = uq_pcset(nord, ndim, pc_type);

% get the NISP quadrature (full tensor)
no1DNodes = 6;
quadrature = uq_quadrature(pcdata.ndim, no1DNodes, pcdata.pc_type);

% get the NISP matrix
K = uq_getNISP(pcdata, quadrature);

% compute system trajectory at quad nodes
[Urlz Vrlz Wrlz] = runNISPRealizations(quadrature);

% save the NISP realizations
save('NISP_data', 'Urlz', 'Vrlz', 'Wrlz');

% time mesh
Tf = 1000;
T = 0 : 1e-2 : Tf;
Nt = length(T);


% get the PC coeffs
U = K * Urlz;
V = K * Vrlz;
W = K * Wrlz; 

% transpose (just for consistency in matrix operations)
U = U';
V = V';
W = W';

% save the results
save('SolutionNISP', 'U', 'V', 'W', 'T', 'pcdata')

% postprocessing 
postprocess(U, V, W, T, pcdata);

