function [xi f] = plotPDF(pcdata, T, tym, X)
% plots the PDF of a given PCE representation 
% at the specified time in the time series.
% The input argument X is a matrix that 
% contains the time evolution of modes 
% of a given system species

[unused itym] = min(abs(T - tym));

pce = X(itym, :); 

Nsamp = 1e5;
U = uq_sample(pcdata, pce, Nsamp);

[f xi] = ksdensity(U);

plot(xi, f, 'linewidth', 2);
set(gca, 'fontsize', 18);
title(['PDF at t = ' num2str(tym)]);

