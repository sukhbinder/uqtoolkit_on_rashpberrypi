function postprocess(U, V, W, T, pcdata)
% Basic postprocessing routine that shows several 
% solution properties. 

% plot the mean and variance of the variabes over time
figure(1);
hold on
plotmeanvar(pcdata, T, U, 'r');
plotmeanvar(pcdata, T, V, 'g');
plotmeanvar(pcdata, T, W, 'b');
set(gca, 'fontsize', 20);
xlabel('Time [-]');
ylabel('Species Mass Fraction [-]');
legend('u', 'v', 'w');

% plot the modes of u over time
figure(2);
plot(T, U, 'linewidth', 2);
set(gca, 'fontsize', 18);
legend({'u_0', 'u_1', 'u_2', 'u_3', 'u_4', 'u_5'}, 'orientation', 'horizonthal', ...
                                                   'location', 'northoutside');
xlim([0 1000])

% plot PDF of u at t = 803
figure(3)
tym = 803;
plotPDF(pcdata, T, tym, U);

% ----------------- plotting subfunction
function plotmeanvar(pcdata, T, X, colstr)

[mu sigma2] = uq_meanvar(pcdata, X);
stddev = sqrt(sigma2);
I = 1 : 1000 : 100000;
errorbar(T(I), mu(I), stddev(I), colstr, 'linewidth', 2);
xlim([0 1000])

