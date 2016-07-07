function [T U V W pcdata] = driver 
% define the path to MATLAB UQ ToolBox 
% modify as needed.
addpath('../../../src_matlab/'); 

% call the Galerkin solver to get the PC modes of the solution
Tf = 1000;
[U V W T pcdata] = solve_dae_ISP(Tf);



% save the results
save('SolutionNISP', 'U', 'V', 'W', 'T', 'pcdata')

% postprocessing 
postprocess(U, V, W, T, pcdata);

