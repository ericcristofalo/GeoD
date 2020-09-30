%% Measurements

if ( s.gt_bool && ~strcmp(init.pose,'experiment') )
    % generate noisy measurements from ground truth data
    y.kappa = deg2rad(y.ang_error); % rotation distribution concentration parameter
    y.y_tau = y.tau; % actual noise parameter for generating measurements
    y.y_kappa = deg2rad(y.ang_error); % actual noise parameter for generating measurements
    y.Tau = zeros(s.n,s.n); % matrix of tau's
    y.Kappa = zeros(s.n,s.n); % matrix of kappa's
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        ind_i = (1:s.d)+(i-1)*s.d;
        ind_j = (1:s.d)+(j-1)*s.d;
        y_t_gt = rob(i).R(:,:,tInd)'*( rob(j).t(:,tInd) - rob(i).t(:,tInd) ); % ground truth
        t_noise = normrnd(0.0,y.y_tau,[3,1]);
        rob(i).y_t(:,j,tInd) = y_t_gt + t_noise; % additive Gaussian noise to translation
        y_R_gt = rob(i).R(:,:,tInd)'*rob(j).R(:,:,tInd); % ground truth
        R_noise = expMap(y.y_kappa*randn(s.d,1));
        rob(i).y_R(:,ind_j,tInd) = R_noise*y_R_gt;
        y.Kappa(i,j) = y.kappa;
        y.Tau(i,j) = y.tau;
    end
    % measurements for se-sync
    measurements = struct;
    measurements.edges = s.M(:,1:2); % each edge in network
    measurements.R = cell(1,s.m);
    measurements.t = cell(1,s.m);
    measurements.kappa = cell(1,s.m);
    measurements.tau = cell(1,s.m);
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        ind_j = (1:s.d)+(j-1)*s.d;
        measurements.R{1,e} = rob(i).y_R(:,ind_j,tInd);
        measurements.t{1,e} = rob(i).y_t(:,j,tInd);
        measurements.kappa{1,e} = y.kappa;
        measurements.tau{1,e} = y.tau;
    end
    
else
    
    % save measurements from data
    y.Tau = zeros(s.n,s.n); % matrix of tau's
    y.Kappa = zeros(s.n,s.n); % matrix of kappa's
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        ind_j = (1:s.d)+(j-1)*s.d;
        rob(i).y_R(:,ind_j,tInd) = measurements.R{1, e};
        rob(i).y_t(:,j,tInd) = measurements.t{1, e};
        y.Kappa(i,j) = measurements.kappa{1, e};
        y.Tau(i,j) = measurements.tau{1, e};
    end
    
    % generate undirected graph from measurements
    measurements_orig = measurements;
    new_edges = measurements.edges;
    for e = 1:s.m
        i = measurements.edges(e, 1);
        j = measurements.edges(e, 2);
        % check if undirected edge exists
        ind_j_check = find(measurements.edges(:, 1)==j);
        ind_i_check = measurements.edges(ind_j_check, 2);
        if ~any(ind_i_check==i)
            % add opposite edge direction for SE-Sync
            disp(['Adding edge: ',num2str(i),' and ',num2str(j)]);
            i_new = size(measurements.edges,1) + 1;
            measurements.edges = [measurements.edges; [j, i]];
            Rji = measurements.R{1, e}';
            measurements.R{1, end+1} = Rji;
            measurements.t{1, end+1} = -Rji*measurements.t{1, e};
            measurements.kappa{1, end+1} = measurements.kappa{1, e};
            measurements.tau{1, end+1} = measurements.tau{1, e};
            % add opposite edge direction for distributed estimation
            new_edges = [new_edges; [j, i]];
            ind_i = (1:s.d)+(i-1)*s.d;
            rob(j).y_R(:,ind_i,tInd) = Rji;
            rob(j).y_t(:,i,tInd) = -Rji*measurements.t{1, e};
            y.Kappa(j,i) = measurements.kappa{1, e};
            y.Tau(j,i) = measurements.tau{1, e};
            s.m = s.m+1;
            
        end
    end
    %    s.m = size(measurements.edges, 1);
    % communication graph properties from measurement data
    s.M = zeros(s.m,3);
    s.M(:,1:2) = new_edges;
    g.L = zeros(s.n,s.n);
    for e = 1:s.m
        i = new_edges(e,1);
        j = new_edges(e,2);
        g.L(i,j) = -1;
        g.L(i,i) = g.L(i,i) + 1;
        s.M(e,3) = sub2ind([s.n,s.n], i, j);
    end
    s.A = (g.L~=0)-eye(s.n);
    eigVal = sort(eig(g.L),'ascend');
    s.fiedler = eigVal(2,1);
    
end


%% Average Measurements to Account for Noise Disagreement
% Applied only to data for distributed consensus estimation. 
% The raw measurements for evaluating the cost are NOT modified. 

measurements.R_pair = cell(1,s.m);
for e = 1:s.m
    i = s.M(e,1);
    j = s.M(e,2);
    ind_i = (1:s.d)+(i-1)*s.d;
    ind_j = (1:s.d)+(j-1)*s.d;
    % Rotation
    yRij = rob(i).y_R(:,ind_j,tInd);
    yRji = rob(j).y_R(:,ind_i,tInd);
    Rd = yRij'*yRji';
    xd = logMap(Rd)/2;
    Rd_ = expMap(xd);
    yRij = yRij*Rd_;
    yRji = Rd_*yRji;
    % Translation
    yt_ij = 0.5*...
        (rob(i).y_t(:,j,tInd) - yRij*rob(j).y_t(:,i,tInd));
    yt_ji = 0.5*...
        (rob(j).y_t(:,i,tInd) - yRji*rob(i).y_t(:,j,tInd));
    % save measurements
    rob(i).y_R(:,ind_j,tInd) = yRij;
    rob(j).y_R(:,ind_i,tInd) = yRji;
%     rob(i).y_t(:,j,tInd) = yt_ij; % leave translation measurements
%     rob(j).y_t(:,i,tInd) = yt_ji;
    measurements.R_pair{1,e} = yRij;

end
disp('Pairwise consistency achieved. ');

