% file: driver.m
%
% purpose: 
%    This is a driver script for galerkin_solver;
%    calls galerkin_solver.m to solve the ODE
%    and then performs some postprocessing. 


[t,Y,pcdata] = galerkin_solver;

% compute the first and second moments:
[mu sigma2] = uq_meanvar(pcdata, Y);

% plot the mean +/- 2 std deviations
figure(1)
stddev = sqrt(sigma2);
%plot(t, mu, 'linewidth', 2);
%hold on;
%plot(t, mu + 2 * stddev, '--r', t, mu - 2 * stddev, '--r');

f = [mu+2*stddev; flipdim(mu-2*stddev,1)];
fill([t(:); flipdim(t(:),1)], f, 'cyan', 'EdgeColor', 'r');
hold on
plot(t, mu, '-k', 'linewidth', 2);


pdfplot = 0;
if pdfplot == 1
   % -----------------------------------------
   % compute and plot the distribution of the 
   % trajectory at t = 1;
   % -----------------------------------------
   figure(2)
   [unused idx] = min(abs(t - 1));   % find the index of time-step 
   pce = Y(idx, :);
   disp('Sampling, please wait ...');
   U = uq_sample(pcdata, pce, 1e4);


   [f,xi] = ksdensity(U);
   plot(xi, f, 'linewidth', 2);
end
