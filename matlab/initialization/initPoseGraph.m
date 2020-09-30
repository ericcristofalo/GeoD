init.poseMean = zeros(3,1);

% load pre-built se-sync data set
if strcmp(init.pose,'sesync_dataset')
    % 3d datasets include: 
    % {'sphere2500','torus3D','grid3D','parking-garage','cubicle','rim'}
    % 2d datasets include: 
    % {'CSAIL', 'manhattan', 'city10000', 'intel', 'ais2klinik', 
    %   'kitti_02'}
    init.dataset = 'parking-garage';
    init.data_dir = [sesync_path,'/data/'];
    init.data_dir = strcat(init.data_dir, init.dataset, '.g2o');
    measurements = load_g2o_data(init.data_dir);
    s.n = max(max(measurements.edges));
    s.m = size(measurements.edges,1);
    % initialize rob data structure
    s.d = length(measurements.t{1});
    s.bound = 0; % desired boundaries of the environment
    s.gt_bool = false; % is ground truth available?
    init.rComm = 1.0; % desired communication radius
    rob = initRobDataStructure(s);

% load T-RO experimental data
elseif strcmp(init.pose,'experiment')
    mydir = pwd;
    idcs = strfind(mydir,'/');
    mydir = mydir(1:idcs(end)-1);
    experiment_data_path = [mydir,'/data/experiment/'];
  
    % read raw camera data
    % poses of robot body frames represented in the global frame
    camera_info = readtable([experiment_data_path,'camera_info.csv']);
    s.n = length(camera_info.globalIndex);
    s.gt_bool = 1; % is ground truth available?
    s.d = 3;
    R_rc = [0,0,1.0; -1.0,0,0; 0,-1.0,0];
    t_r_c = [0.05, 0.0, 0.07]';
    rob = initRobDataStructure(s);
    for i = 1:s.n
        rob(i).t = ...
            [camera_info.tx(i); camera_info.ty(i); camera_info.tz(i)];
        R = [camera_info.r11(i), camera_info.r12(i), camera_info.r13(i);
            camera_info.r21(i), camera_info.r22(i), camera_info.r23(i);
            camera_info.r31(i), camera_info.r32(i), camera_info.r33(i)];
        [U,S,V] = svd(R);
        rob(i).R = U*V'; % constrain data to rotation matrix
        rob(i).image_name = camera_info.imagePath(i);
    end
    
    % read raw measurement data (in local camera frame)
    measurement_info = readtable([experiment_data_path,...
        'measurement_info.csv']);
    s.bound = 0; % desired boundaries of the environment
    init.rComm = 0.75; % desired communication radius
    y.Tau = zeros(s.n,s.n); % matrix of tau's
    y.Kappa = zeros(s.n,s.n); % matrix of kappa's
    s.m = length(measurement_info.globalIndexI);
    measurements = struct;
    measurements.edges = ...
        [measurement_info.globalIndexI + 1, ...
         measurement_info.globalIndexJ + 1];
    measurements.R = cell(1,s.m);
    measurements.t = cell(1,s.m);
    measurements.tau = cell(1,s.m);
    measurements.kappa = cell(1,s.m);
    max_num_features = max(measurement_info.numberOfMatches);
    for e = 1:s.m
        i = measurements.edges(e,1);
        j = measurements.edges(e,2);
        
        % read raw relative camera measurements
        % relative poses of camera body frame i represented camera body frame j
        t_cj_ci = ...
            [measurement_info.tx(e); measurement_info.ty(e); measurement_info.tz(e)];
        R_cjci = ...
            [measurement_info.r11(e), measurement_info.r12(e), measurement_info.r13(e);
            measurement_info.r21(e), measurement_info.r22(e), measurement_info.r23(e);
            measurement_info.r31(e), measurement_info.r32(e), measurement_info.r33(e)];
        [U,S,V] = svd(R_cjci); % constrain data to rotation matrix
        R_cjci = U*V';
        
        % compare to ground truth
        R_ri = rob(i).R;
        R_rj = rob(j).R;
        R_rirj_gt = rob(i).R'*rob(j).R;
        R_cicj_gt = R_rc'*R_rirj_gt*R_rc;
        R_cjci_gt = R_cicj_gt';
        t_ri = rob(i).t;
        t_rj = rob(j).t;
        t_ci_gt = t_ri + rob(i).R*t_r_c;
        t_cj_gt = t_rj + rob(j).R*t_r_c;
        t_ci_cj_gt = R_rc'*rob(i).R'*(t_cj_gt - t_ci_gt);
        t_cj_ci_gt = R_rc'*rob(j).R'*(t_ci_gt - t_cj_gt);
        
        % convert ji measurement to ij for my convention
        R_cicj = R_cjci';
        R_rirj = R_rc*R_cicj*R_rc';
        t_ci_cj = -R_cicj_gt*t_cj_ci;
        t_ri_rj = R_rc*t_ci_cj + R_rirj_gt*t_r_c + t_r_c;

        measurements.t{e} = t_ri_rj;
        measurements.R{e} = R_rirj;
        meas_uncert = 1.0/(measurement_info.numberOfMatches(e)/max_num_features);
        measurements.tau{1,e} = meas_uncert;
        measurements.kappa{1,e} = meas_uncert;
        y.Tau(i,j) = meas_uncert;
        y.Kappa(i,j) = meas_uncert;
    end
    s.M = [measurements.edges, zeros(s.m,1)];
    clear measurement_info;

% generate random initial poses
elseif ( strcmp(init.pose,'line') || strcmp(init.pose,'random') )
    % initialize rob data structure
    s.d = 3; % dimension of the problem
    s.n = 25; % number of robots
    s.bound = 5; % desired boundaries of the environment
    s.gt_bool = true; % is ground truth available?
    init.rComm = 1.0; % desired communication radius
    rob = initRobDataStructure(s);
    for i = 1:s.n
        % generated poses
        if strcmp(init.pose,'line') % line initial positions
            rob(i).t(:,tInd) = 2.0*s.bound*[rand(1,1);0;0] - s.bound;
        elseif strcmp(init.pose,'random') % random initial positions
            rob(i).t(:,tInd) = 2.0*s.bound*rand(3,1) - s.bound;
        end
        % generated rotations
        rob(i).R(:,:,tInd) = randRot(s.d);
        % ensure communication with at least one neighbor
        if ( i>1 )
            neighInd = ceil((i-1)*rand(1,1));
            d = 0.0;
            rad = (init.rComm-d)*rand(1,1)+d; % uniform random radius from previous robot
            relPose = rob(i).t(:,tInd) - rob(neighInd).t(:,tInd);
            rob(i).t(:,tInd) = rob(neighInd).t(:,tInd) + rad*relPose/norm(relPose);
        end
        % count mean position
        init.poseMean = init.poseMean + rob(i).t(:,tInd); % mean position
    end

% generate circle formation
elseif strcmp(init.pose,'circle')
    % initialize rob data structure
    s.d = 3; % dimension of the problem
    s.n = 25; % number of robots
    s.bound = 5; % desired boundaries of the environment
    s.gt_bool = true; % is ground truth available?
    init.rComm = 4.0; % desired communication radius
    init.radius = 2.6;
    alpha = 2*pi/s.n;
    rob = initRobDataStructure(s);
    for i = 1:s.n
        % generated poses
        Rad = euler2rot([0;0;(i-1)*alpha]);
        rob(i).t(:,tInd) = Rad*init.radius*[1;0;0];
        rob(i).R(:,:,tInd) = randRot(s.d);
        % count mean position
        init.poseMean = init.poseMean + rob(i).t(:,tInd); % mean position
    end

% generate grid formation
elseif strcmp(init.pose,'grid')
    % initialize rob data structure
    s.d = 3; % dimension of the problem
    init.numRobotsPerSide = 3; % number of robots per side
    s.n = (init.numRobotsPerSide)^s.d; % number of robots
    s.gt_bool = true; % is ground truth available?
    init.rComm = 1.66; % desired communication radius
    s.bound = init.rComm*init.numRobotsPerSide*1.2; % desired boundaries of the environment
    init.side = round((s.n)^(1/3));
    init.sqSize = init.rComm*init.side/1;
    rob = initRobDataStructure(s);
    for i = 1:init.side
        for ii = 1:init.side
            for iii = 1:init.side
                % position
                ind = (init.side*init.side*(i-1) + init.side*(ii-1) + iii);
                rob(ind).t(:,tInd) = [init.sqSize/init.side*i,...
                    init.sqSize/init.side*ii,...
                    init.sqSize/init.side*iii];
                % rotation
                rob(ind).R(:,:,tInd) = randRot(s.d);
                % count mean position
                init.poseMean = init.poseMean + rob(ind).t(:,1);
            end
        end
    end

% generate sphere formation
% https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
elseif strcmp(init.pose,'sphere')
    % initialize rob data structure
    s.d = 3; % dimension of the problem
    s.n = 50; % number of robots
    s.bound = 5; % desired boundaries of the environment
    s.gt_bool = true; % is ground truth available?
    init.radius = 3.0;
    init.rComm = 3.0;
    init.a = 4.0*pi*(1.0)^2/s.n;
    init.d = sqrt(init.a);
    init.m_alpha = round(pi/init.d);
    init.d_alpha = pi/init.m_alpha;
    init.d_beta = init.a/init.d_alpha;
    ind = 0;
    rob = initRobDataStructure(s);
    for i = 0:(init.m_alpha-1)
        init.alpha = pi*(i+0.5)/init.m_alpha;
        init.m_beta = round(2*pi*sin(init.alpha)/init.d_beta);
        for ii = 0:(init.m_beta-1)
            % current robot index
            ind = ind+1;
            % check final robot count
            if (ind>s.n)
                break;
            end
            init.beta = 2*pi*ii/init.m_beta;
            % generated poses
            rob(ind).t(:,tInd) = init.radius*...
                [sin(init.alpha)*cos(init.beta);
                sin(init.alpha)*sin(init.beta);
                cos(init.alpha)];
            % generated rotations
            rob(ind).R(:,:,tInd) = randRot(s.d);
            % count mean position
            init.poseMean = init.poseMean + rob(ind).t(:,tInd); % mean position
        end
    end
    % duplicate extra robots if required
    if (ind<s.n)
        temp = s.n-ind;
        for i = 1:temp
            rob(ind+i).t(:,tInd) = rob(ind).t(:,tInd);
            rob(ind+i).R(:,:,tInd) = rob(ind).R(:,:,tInd);
        end
    end

else
    
    error('Error in initPoseGraph.m: incorrect pose graph selection');
    
end


%% Finalize Initialization of Graph Structures

if ( s.gt_bool && ~strcmp(init.pose,'experiment') )
    % center initial conditions
    init.poseMean = init.poseMean/s.n;
    for i = 1:s.n
        rob(i).t(:,tInd) = rob(i).t(:,tInd) - init.poseMean + s.bound/2.0*ones(3,1);
    end
    % communication graph initialization
    A = zeros(s.n,s.n);
    for i = 1:s.n
        for j = i:s.n
            if i~=j
                testRad = norm( rob(i).t(:,tInd) - rob(j).t(:,tInd) );
                if testRad<=init.rComm
                    A(i,j) = 1; % unweighted graph edges
                end
            end
        end
    end
    s.A = A+A'; % adjacency matrix
    D = diag(sum(s.A,2)); % degree matrix
    g.L = D-s.A; % Laplacian matrix
    [ind_i,ind_j] = find(s.A); % components of directed connections
    s.M = [ind_i,ind_j,find(s.A)]; % matrix of edges
    % [i index, j index, index of L matrix];
    s.m = size(s.M,1);
    eigVal = sort(eig(g.L),'ascend');
    s.fiedler = eigVal(2,1);
    
else
    
    % communication graph extraction from data
    s.M = zeros(s.m,3);
    s.M(:,1:2) = measurements.edges;
    
    % extract graph laplacian from data
    g.L = zeros(s.n,s.n);
    for e = 1:s.m
        i = measurements.edges(e,1);
        j = measurements.edges(e,2);
        g.L(i,j) = -1;
        g.L(i,i) = g.L(i,i) + 1;
        s.M(e,3) = sub2ind([s.n,s.n], i, j);
    end
    s.A = (g.L~=0)-eye(s.n);
    eigVal = sort(eig(g.L),'ascend');
    s.fiedler = eigVal(2,1);
    
end
