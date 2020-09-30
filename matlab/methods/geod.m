%% Initialization

tInd = 1;

% initialize distributed estimation data structures
if ~isfield(rob, 't_hat_d')
    for i = 1:s.n
        rob(i).t_hat_d = zeros(s.d,s.steps);
        rob(i).R_hat_d = repmat(eye(s.d),1,1,s.steps);
    end
end
s.mle_error_d = zeros(3,s.steps);
s.dist_time = zeros(1,s.steps);

% initialize estimates
if ( ~isfield(init,'R') )
    [t_hat,R_hat] = initPoseEst(s.initType,rob,s,measurements,y);
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        rob(i).t_hat_d(:,tInd) = t_hat(:,i);
        rob(i).R_hat_d(:,:,tInd) = R_hat(:,ind_i);
        init.t(:,i) = t_hat(:,i);  % save for other comparisons
        init.R(:,ind_i) = R_hat(:,ind_i);  % save for other comparisons
    end
else
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        rob(i).t_hat_d(:,tInd) = init.t(:,i);
        rob(i).R_hat_d(:,:,tInd) = init.R(:,ind_i);
    end
end


%% GeoD Simulation

disp(' ');
disp('========== GEOD ==========');
disp(' ');

for tInd = 2:s.steps
    
    disp(['--- Time step ',num2str(tInd),' of ',num2str(s.steps),' ---']);
    
    % update robot poses here if the pose graph is dynamically changing
    % (static robot poses in this version)
    
    % distributed estimation control step
    tic
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        ind_m = find(s.M(:,1)==i);
        m = size(ind_m,1);
        ind_neighs = s.M(ind_m,2);
        gamma = s.dt;
        if ( m~=0 && tInd~=1 ) % if connected to neighbors
            
            % ---------- rotation pose graph consensus ----------
            if strcmp(distOption,'geod')
                
                % geod rotation update (equation 3b)
                omega = zeros(s.d,1);
                for e = ind_m'
                    j = s.M(e,2);
                    ind_j = (1:s.d)+(j-1)*s.d;
                    omega = omega + ...
                        logMap(rob(i).R_hat_d(:,:,tInd-1)'* ...
                               rob(j).R_hat_d(:,:,tInd-1)* ...
                               rob(i).y_R(:,ind_j,1)');
                end
                % rotation dynamics update
                rob(i).R_hat_d(:,:,tInd) = ...
                    rob(i).R_hat_d(:,:,tInd-1)*expMap(gamma*omega);

            elseif strcmp(distOption,'angle_axis')
                
                % angle-axis update (2019 cdc paper)
                omega = zeros(s.d,1);
                for e = ind_m'
                    j = s.M(e,2);
                    ind_j = (1:s.d)+(j-1)*s.d;
                    omega = omega + ...
                        logMap(rob(j).R_hat_d(:,:,tInd-1)) - ...
                        logMap(rob(i).R_hat_d(:,:,tInd-1)) - ...
                        logMap(rob(i).y_R(:,ind_j,1));
                end
                % angle-axis dynamics
                omega = gamma*omega;
                x_t = logMap(rob(i).R_hat_d(:,:,tInd-1));
                theta = norm(x_t);
                beta = 0;
                x_t_ = zeros(s.d,1);
                if (theta~=0)
                    x_t_ = x_t./theta;
                    beta = 1 - 0.25*(theta*sin(theta)/(sin(0.5*theta)^2));
                end
                x_exp = hat(x_t_);
                L_x_i = eye(s.d) + 0.5*theta*x_exp + beta*x_exp*x_exp;
                % angle-axis discrete differentiation
                axis_new = x_t + s.dt*(L_x_i*omega);
                rob(i).R_hat_d(:,:,tInd) = expMap(axis_new);
                
            else
                
                % perform no rotation update
                rob(i).R_hat_d(:,:,tInd) = ...
                    rob(i).R_hat_d(:,:,tInd-1);
                if ( tInd==2 && i==1 )
                    disp('Warning in geod.m: no rotation update performed')
                end
                
            end
            % ---------- ---------- ---------- ---------- ----------
            
            % ---------- translation pose graph consensus ----------
            % geod translation update (equation 3a)
            nu = zeros(s.d,1); % translation
            for e = ind_m'
                j = s.M(e,2);
                % translation consensus term in world frame updates
                % with built-in measurement averaging
                nu = nu + ...
                    rob(j).t_hat_d(:,tInd-1) - rob(i).t_hat_d(:,tInd-1) ...
                    - 0.5*rob(i).R_hat_d(:,:,tInd-1)*rob(i).y_t(:,j,1) ...
                    + 0.5*rob(j).R_hat_d(:,:,tInd-1)*rob(j).y_t(:,i,1);
            end
            % translation dynamics
            rob(i).t_hat_d(:,tInd) = rob(i).t_hat_d(:,tInd-1) + gamma*nu;
            % ---------- ---------- ---------- ---------- ----------
            
        end % if connected to neighbors
    end % for each robot
    s.dist_time(1,tInd) = toc;

end
disp(' ');
disp('===== END GEOD =====');
disp(' ');


%% Compute Distributed Pose Graph Error

for tInd_ = 1:tInd
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
            norm(rob(j).t_hat_d(:,tInd_) - rob(i).t_hat_d(:,tInd_) - ...
                 rob(i).R_hat_d(:,:,tInd_)*measurements.t{e},2).^2;
        R_error_sesync = ...
            norm(rob(j).R_hat_d(:,:,tInd_) - ...
                 rob(i).R_hat_d(:,:,tInd_)*measurements.R{e},'fro').^2;
        R_error_geodesic = norm(logMap( ...
            rob(i).R_hat_d(:,:,tInd_)'*rob(j).R_hat_d(:,:,tInd_)* ...
            measurements.R{e}'),2).^2;
        s.mle_error_d(1,tInd_) = s.mle_error_d(1,tInd_) + ...
                                 measurements.tau{1,e}*t_error;
        s.mle_error_d(2,tInd_) = s.mle_error_d(2,tInd_) + ...
                                 measurements.kappa{1,e}*R_error_sesync;
        s.mle_error_d(3,tInd_) = s.mle_error_d(3,tInd_) + ...
                                 measurements.kappa{1,e}*R_error_geodesic;
    end
    
end

