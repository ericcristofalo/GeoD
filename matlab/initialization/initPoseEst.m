function [t_hat, R_hat] = initPoseEst(type,rob,s,measurements,y)

% centralized structure
t_hat = zeros(s.d,s.n);
R_hat = zeros(s.d,s.n*s.d);

if strcmp(type,'random')
    
    % basic random initialization
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        t_hat(:,i) = normrnd(0.0,1.0,[s.d,1]);
        R_hat(:,ind_i) = randRot;
    end
    
    % ensure distributed estimation condition is met
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        ind_m = find(s.M(:,1)==i);
        for e = ind_m'
            j = s.M(e,2); % neighbor index
            ind_j = (1:s.d)+(j-1)*s.d;
            init_check = logMap(R_hat(:,ind_i)'*R_hat(:,ind_j)*measurements.R{e}');
            while norm(init_check, 2) < (pi - 0.1)
                R_hat(:,ind_i) = randRot;
                init_check = logMap(R_hat(:,ind_i)'*R_hat(:,ind_j)*measurements.R{e}');
            end
        end
    end
    
elseif strcmp(type, 'identity')
    
    % Basic Initialization
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        t_hat(:,i) = zeros(s.d, 1);
        R_hat(:,ind_i) = eye(s.d);
    end
    
elseif strcmp(type,'gps')
    
    % initialize estimates to noisy gps measurements
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        % measurement dependent noise
        t_noise = normrnd(0.0,y.tau,[s.d,1]);
        R_noise = expMap(normrnd(0.0,deg2rad(y.ang_error),[s.d,1]));
        %       % fixed gps noise
        %       t_noise = normrnd(0.0,3.0,[s.d,1]);
        %       R_noise = expMap(normrnd(0.0,deg2rad(30.0),[s.d,1]));
        % generate noisy measurements
        t_hat(:,i) = rob(i).t(:,1) + t_noise;
        R_hat(:,ind_i) = rob(i).R(:,:,1)*R_noise;
    end
    
elseif strcmp(type,'chordal')
    
    % initialize translation estimates from spanning tree through network
    rootInd = 1;
    s = buildSpanningTree(s,rootInd);
    [t_hat,~] = spanningTreeTransformations(s,rob,rootInd,1,1);
%     % basic random initialization
%     for i = 1:s.n
% %         t_hat(:,i) = normrnd(0.0,10.0,[s.d,1]);
%         t_hat(:,i) = (2*rand(s.d,1)-1)*s.bound;
%     end
%     % generate noisy gps measurements
%     for i = 1:s.n
%         t_noise = normrnd(0.0,y.tau,[s.d,1]);
%         t_hat(:,i) = rob(i).t(:,1) + t_noise;
%     end
    % centralized structure
    R_hat = chordal_initialization(measurements);
    
elseif strcmp(type,'spanning_tree')
    
    % Initialize Estimates from Spanning Tree Through Network
    rootInd = 1;
    s = buildSpanningTree(s,rootInd);
    [t_hat,R_hat] = spanningTreeTransformations(s,rob,rootInd,1,1);
    
elseif strcmp(type,'chordal_distributed')
    
    % Jacobi solution to chordal relaxation
    % convert my structures
    N = s.n;
    d = s.d;
    w = randn(d, N);
    R = zeros(d, d, N);
    for i = 1:N
        R(:, :, i) = expMap(hat(w(:, i)));
    end
    neighbors = cell(N, 1);
    for i = 1:N
        neighbors{i} = find(s.A(i, :));
    end
    Rij = cell(N, N);
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        Rij{i,j} = measurements.R{e};
    end
    diameter = 100; % fixed number of iterations
    T = diameter + 0; % number of iterations
    blk3 = @(W) blkdiag(W, W, W);
    id = @(r) (r-1)*9 + (1:9);
    Hgs = cell(N, N); % H matrix for Gauss-Seidel iterations
    for i = 1:N
        Hgs{i, i} = zeros(d^2);
        for j = neighbors{i}
            Hgs{i, i} = Hgs{i, i} + blk3(Rij{i, j}*Rij{i, j}' + eye(d));
            Hgs{i, j} = blk3(-Rij{i, j} - Rij{j, i}');
        end
    end
    Hcent_cell = Hgs;
    for i = 1:N
        for j = 1:N
            if isempty(Hcent_cell{i, j})
                Hcent_cell{i, j} = zeros(d^2);
            end
        end
    end
    % run
    Rgs = cell(N, T);
    ygs = nan(d^2*N, T+1);
    ygs(:, 1) = [reshape(R(:, :, 1)', [d^2, 1]); zeros(d^2*(N-1), 1)];
    gamm = 1;
    successive = false; % successive = true means updates are asynchronous
    for t = 1:T
        [~, order] = sort(rand(N, 1)); % random order of nodes (only matters for successive version)
        for i = order'
            gi = zeros(9, 1);
            for j = neighbors{i}
                if any(isnan(ygs(id(j), t+1))) || ~successive
                    gi = gi + Hgs{i, j}*ygs(id(j), t);
                else
                    gi = gi + Hgs{i, j}*ygs(id(j), t+1);
                end
            end
            ygs(id(i), t+1) = (1-gamm)*ygs(id(i), t) + gamm*Hgs{i, i}\-gi;
        end
        % projection
        for i = 1:N
            [Utemp, ~, Vtemp] = svd(reshape(ygs(id(i), t), 3, 3)');
            Rgs{i, t} = Utemp*diag([1 1 det(Utemp*Vtemp')])*Vtemp';
        end  
    end
    % save data
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        R_hat(:,ind_i) = Rgs{i,T};
    end
    % basic random initialization for translation
    for i = 1:s.n
%         t_hat(:,i) = normrnd(0.0,s.bound/2,[s.d,1]);
%         t_hat(:,i) = (2*rand(s.d,1)-1)*s.bound/2;
        t_hat(:,i) = zeros(s.d,1);
    end
%     % generate noisy gps measurements
%     for i = 1:s.n
%         t_noise = normrnd(0.0,y.tau,[s.d,1]);
%         t_hat(:,i) = rob(i).t(:,1) + t_noise;
%     end
    
else
    
    error('Error in initPoseEst.m: incorrect type of initialization provided');

end

end