function plottraj(T, X)
% used to make a plot of the system trajectory

plot(T, X(:,1), T, X(:,2), T, X(:,3), 'linewidth', 2); 
xlim([-5 1000]); 
ylim([0 1]); 
set(gca, 'fontsize', 20);
xlabel('Time [-]');
ylabel('Species Mass Fraction [-]');
legend('u', 'v', 'w');
