%% SE-Sync

tic;

% set manopt options (if desired)
Manopt_opts.tolgradnorm = 1e-10;  % Stopping tolerance for norm of Riemannian gradient
Manopt_opts.rel_func_tol = 1e-10;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value
Manopt_opts.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxiter = 100;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxinner = 100;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
%manopt_options.maxtime = 60*60;  % Maximum computation time to allow, in seconds
%manopt_options.solver = @steepestdescent;  % Select Manopt solver to use: {trustregions (default), conjugategradient, steepestdescent}
Manopt_opts.verbosity = 0;

% set se-sync options (if desired)
SE_Sync_opts.r0 = 3;  % Initial maximum-rank parameter at which to start the Riemannian Staircase
SE_Sync_opts.rmax = 10;  % Maximum maximum-rank parameter at which to terminate the Riemannian Staircase
SE_Sync_opts.eig_comp_rel_tol = 1e-10;  % Relative tolerance for the minimum-eigenvalue computation used to test for second-order optimality with MATLAB's eigs() function
SE_Sync_opts.min_eig_lower_bound = -1e-10;  % Minimum eigenvalue threshold for accepting a maxtrix as numerically positive-semidefinite
SE_Sync_opts.Cholesky = false;  % Select whether to use Cholesky or QR decomposition to compute orthogonal projections
use_chordal_initialization = true;  % Select whether to use the chordal initialization, or a random starting point

% initialization
if exist('measurements_orig', 'var')
    [t_hat,R_hat] = initPoseEst(s.initType,rob,s,measurements_orig,y);
else
    [t_hat,R_hat] = initPoseEst(s.initType,rob,s,measurements,y);
end
init.t = t_hat;
init.R = R_hat;
Y0 = vertcat(R_hat, zeros(SE_Sync_opts.r0 - s.d, s.n*s.d));

% run se-sync
if exist('measurements_orig', 'var')
    [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = ...
        SE_Sync(measurements_orig, Manopt_opts, SE_Sync_opts, Y0);
else
    [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = ...
        SE_Sync(measurements, Manopt_opts, SE_Sync_opts, Y0);
end

% save results
for i = 1:s.n
    ind_i = (1:s.d)+(i-1)*s.d;
    rob(i).t_hat = xhat.t(:,i);
    rob(i).R_hat = xhat.R(:,ind_i);
end
s.Fxhat = Fxhat; % objective value
s.SDPval = SDPval; % sdp optimal value
s.cent_time = toc;


%% Compute Centralized Pose Graph Error

s.mle_error = zeros(3,1); % mle error structure
if exist('measurements_orig', 'var')
    num_edges = min(size(measurements.edges,1), ...
        size(measurements_orig.edges,1));
else
    num_edges = size(measurements.edges,1);
end
for e = 1:num_edges
    i = measurements.edges(e,1);
    ind_i = (1:s.d)+(i-1)*s.d;
    j = measurements.edges(e,2);
    ind_j = (1:s.d)+(j-1)*s.d;
    t_error = ...
        norm(rob(j).t_hat - rob(i).t_hat - ...
             rob(i).R_hat*measurements.t{e},2).^2;
    R_error_sesync = ...
        norm(rob(j).R_hat - rob(i).R_hat*measurements.R{e},'fro').^2;
    R_error_geodesic = norm(logMap(...
        rob(i).R_hat'*rob(j).R_hat*measurements.R{e}')...
        ,2).^2;
    s.mle_error(1,:) = s.mle_error(1,:) + measurements.tau{1,e}*t_error;
    s.mle_error(2,:) = s.mle_error(2,:) + ...
                       measurements.kappa{1,e}*R_error_sesync;
    s.mle_error(3,:) = s.mle_error(3,:) + ...
                       measurements.kappa{1,e}*R_error_geodesic;
end

