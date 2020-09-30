%% Initialization

clear; clc; close all;

rng(0);
folders_to_add = ...
    {'initialization', 'math', 'methods', 'plotting', 'transformations'};
for i = 1:length(folders_to_add)
    addpath(genpath( folders_to_add{i} ));
end

% simulation time steps
s.dt = 0.025;
s.steps = 101;
s.tEnd = s.dt*s.steps;
s.tspan = linspace(0,s.steps,s.steps)*s.tEnd/s.steps;
tInd = 1;

% set local path to se-sync if using standard slam datasets or se-sync
% algorithm
% SET SE-SYNC PATH HERE:
% sesync_path = '/path/to/SE-Sync';
% addpath(genpath(sesync_path));
sesync_comparison = false;

% initialize pose graph
% pose graph options include:
% {'line', 'random', 'circle', 'grid', 'sphere', 'sesync_dataset',
%	'experiment'};
% The first five graphs are simulated pose graphs that can be modified in
%	initPoseGraph.m. 
% The sesync_data set options chooses a standard SLAM dataset included with
%   SE-Sync. SE-Sync must be downloaded and installed. The path must be set
%   above. 
% The experimental pose graph represents the multi-UAV dataset used in the
%   paper. 
init.pose = 'grid';
initPoseGraph;

% initialization pose graph optimization
% initialization options include:
% {'identity', 'random', 'gps', 'chordal', 'spanning_tree'};
s.initType = 'gps';
if strcmp(init.pose,'sesync_dataset') && strcmp(s.initType,'gps')
   error('Error in main.m: gps initialization can not be used with sesync datasets');
end

% generate measurements if not using dataset
y.tau = 0.5; % translation Gaussian error (meters)
y.ang_error = 30.0; % Langevin covariance for rotation error (degrees)
initMeasurements;


%% SE-Sync

if sesync_comparison
    sesync;
end


%% GeoD

% consensus-based distributed pose graph options include:
% {'geod', 'angle_axis'};
distOption = 'geod';
geod;


%% Display Results

disp('========== FINAL RESULTS ==========');
disp(' ');

% centralized results
if sesync_comparison
    disp('=== SE-Sync ===');
    disp(' ');
    disp('Did SE-Sync reach certified global chordal minimum?')
    if abs(s.SDPval(1,1)-s.Fxhat(1,1))<1E-10
        disp('    Yes');
    else
        disp('    No');
    end
    disp('Total time (seconds): ');
    fprintf('    %.3f\n', SE_Sync_info.total_computation_time + s.cent_time);
    disp('Chordal value: ');
    fprintf('    %.3f\n', s.mle_error(1,1)+s.mle_error(2,1));
    disp('Geodesic value: ');
    fprintf('    %.3f\n', s.mle_error(1,1)+s.mle_error(3,1));
end
disp(' ');

% distributed error results
% compute values at convergence
convergence_threshold = 1E-2; 
mle_total = s.mle_error_d(1,:)+s.mle_error_d(3,:);
mle_total_chord = s.mle_error_d(1,:)+s.mle_error_d(2,:);
conv_ind = abs(diff(mle_total))>convergence_threshold;
mle_total = mle_total(1,conv_ind==1);
disp('=== GeoD ===')
disp(' ');
disp('Total time/agent at convergence (seconds) :');
fprintf('    %.3f\n', sum(s.dist_time(1,1:length(mle_total)))/s.n);
disp('Total GeoD iterations at convergence: ');
fprintf('    %.0f\n', length(mle_total));
disp('Chordal value: ');
fprintf('    %.3f\n', mle_total_chord(length(mle_total)));
disp('Geodesic value: ');
fprintf('    %.3f\n', mle_total(end));

% Plot Estimation from Distributed SE-Sync
plotResults;

