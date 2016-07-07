function [xi f] = plotPDF(pcdata, T, tym, X)
% plots the PDF of a given random variable given its 
% PC representation at the specified time in 
% the time series.  

[unused itym] = min(abs(T - tym));

pce = X(itym, :); 

Nsamp = 1e5;
U = uq_sample(pcdata, pce, Nsamp);

[f xi] = ksdensity(U);

plot(xi, f, 'linewidth', 2);
set(gca, 'fontsize', 18);
title(['PDF at t = ' num2str(tym)]);

