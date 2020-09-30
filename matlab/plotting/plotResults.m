%% Plotting Setup

p.R_plot_w = eye(3); % transformation to world frame
p.cVec = colormap; % default color map
p.cInd = ceil(linspace(1,64,s.n)); % color map scaled to number of robots
p.arrowLength = 0.25; % size of coordinate systems
p.gtAlpha = 1; % translucent value for estimates (0 is transparent), 0.3 normally
p.estAlpha = 1; % translucent value for estimates
p.commColor = [0.5*ones(1,3),p.gtAlpha]; % communication link color, 0.5 for light color
p.commEstColor = [0.5*ones(1,3),p.gtAlpha]; % estimate communication link color
p.maxPlotN = 251;

% final plotting values
p.plotViewValues = [-1,-1.7,0.8];
if strcmp(init.pose, 'random')
%     p.plotViewValues = [71, 30];
elseif strcmp(init.pose, 'circle')
%     p.plotViewValues = [71, 88];
elseif strcmp(init.pose, 'sesync_dataset')
    p.plotViewValues = [-91.7361, 33.1815];
end


%% Plot Initial Ground Truth Environment

if s.gt_bool
    
    figure(1);
    clf(1);
    titleString = 'Ground Truth';
    hold on;
%     % plot world reference frame
%     h = plotCoordSys([0; 0; 0; 0; 0; 0], 'world', 50, [1,1,1], 0, 1, 2);
    
    % plot communication network
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        x_pts = [rob(i).t(1,1); rob(j).t(1,1)];
        y_pts = [rob(i).t(2,1); rob(j).t(2,1)];
        z_pts = [rob(i).t(3,1); rob(j).t(3,1)];
        plot3(x_pts,y_pts,z_pts, ...
            'Color',p.commColor,'LineStyle','-');
    end
    
    % plot robot coordinate systems
    if ( s.n<p.maxPlotN )
        for i = 1:s.n
            % for coordinate system labels
            % plotCoordSys([rob(i).t(:,tInd);rot2euler(rob(i).R(:,:,tInd))], ...
            %              num2str(i), 0, p.cVec(p.cInd(i),:), p.gtAlpha, p.arrowLength, 2);
            % no coordinate system labels
            plotCoordSys([rob(i).t(:,1);rot2euler(rob(i).R(:,:,1))], ...
                '', 1, p.cVec(p.cInd(i),:), p.gtAlpha, p.arrowLength, 2);
        end
    end
    hold off;
    box on;
    rotate3d on;
    axis equal;
    view(p.plotViewValues);
    title(titleString,'Interpreter','LaTex');
    % define plot volume
    p.plotVol = zeros(1,6);
    min_bound = 1E6*ones(3,1);
    max_bound = -1E6*ones(3,1);
    for i = 1:s.n
        ind = find(rob(i).t(:)<min_bound);
        for j = ind'
            p.plotVol((j-1)*2+1) = rob(i).t(j,1);
            min_bound(j) = rob(i).t(j,1);
        end
        ind = find(rob(i).t(:)>max_bound);
        for j = ind'
            p.plotVol((j-1)*2+2) = rob(i).t(j,1);
            max_bound(j) = rob(i).t(j,1);
        end
    end
    offset = 0.1;
    p.plotVol = p.plotVol + [-offset,offset,-offset,offset,-offset,offset];
    axis(p.plotVol);
%     axis off; % to hide axes and background
    
end


%% Plot Centralized Result

if sesync_comparison

    figure(2);
    clf(2);
    titleString = 'SE-Sync Result';
    hold on;

    % plot centralized estimate
    robInd = 1; % plot estimate from robInd's perspective
    i = robInd;
    p.robInd_t = zeros(s.d,1);
    p.robInd_R = eye(s.d);
    if ( s.gt_bool==1 ) % translate estimate if ground truth is available
        p.robInd_t = rob(i).t(:,1);
        p.robInd_R = rob(i).R(:,:,1);
    end

    % robot's estimate to plot
    ind_i = (1:s.d)+(i-1)*s.d;
    ind_m = find(s.M(:,1)==i);
    m = size(ind_m,1);
    p.rob.t_plot = xhat.t; % values directly from se-sync
    p.rob.R_plot = xhat.R; % values directly from se-sync

    % transform estimation about frame i in world frame for comparison
    t_trans = p.robInd_t - p.rob.t_plot(:,i);
    R_trans = p.rob.R_plot(:,ind_i)'*p.robInd_R;
    [p.t_cur, p.R_cur] = ...
        rigidBodyTrans(p.rob.t_plot,p.rob.R_plot,t_trans,R_trans,i);
    % remove agents with no estimate from ith point of view
    zeroInd = find(p.rob.t_plot(1,:)==0); % indices to remove from lack of relative measurement
    zeroInd(zeroInd==i) = [];
    % plot communication network
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        x_pts = [p.t_cur(1,i); p.t_cur(1,j)];
        y_pts = [p.t_cur(2,i); p.t_cur(2,j)];
        if ( s.d==2 )
            plot(x_pts,y_pts, ...
                'Color',p.commEstColor,'LineStyle','-','LineWidth',1);
        else
            z_pts = [p.t_cur(3,i); p.t_cur(3,j)];
            %          plot3(x_pts,y_pts,z_pts, ...
            %             'Color',p.commEstColor,'LineStyle','-','LineWidth',1);
            plot3(x_pts,y_pts,z_pts, ...
                'Color',p.commColor,'LineStyle','-'); % for final figures
        end
    end
    % plot coordinate systems
    if ( s.n<p.maxPlotN )
        for ii = 1:s.n
            if ( any(p.t_cur(:,ii)~=0) && ~any(ii==zeroInd) )
                ind_ii = (1:s.d)+(ii-1)*s.d;
                % for coordinate system labels
                % plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
                %              ['$\hat{',num2str(ii),'}$'], 50, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
                % no coordinate system labels
                plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
                    '', 0, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
            end
        end
    end
    hold off;
    box on;
    rotate3d on;
    axis equal;
    view(p.plotViewValues);
    title(titleString,'Interpreter','LaTex');
    % define plot volume
    p.plotVol = zeros(1,6);
    min_bound = 1E6*ones(3,1);
    max_bound = -1E6*ones(3,1);
    for i = 1:s.n
        ind = find(p.t_cur(:,i)<min_bound);
        for j = ind'
            p.plotVol((j-1)*2+1) = p.t_cur(j,i);
            min_bound(j) = p.t_cur(j,i);
        end
        ind = find(p.t_cur(:,i)>max_bound);
        for j = ind'
            p.plotVol((j-1)*2+2) = p.t_cur(j,i);
            max_bound(j) = p.t_cur(j,i);
        end
    end
    offset = 0.1;
    p.plotVol = p.plotVol + [-offset,offset,-offset,offset,-offset,offset];
    axis(p.plotVol);
    % axis off; % to hide axes and background

end


%% Plot Distributed Result

    
figure(3);
clf(3);
titleString = 'GeoD Result';
hold on;
robInd = 1; % plot estimate from robInd's perspective
p.robInd_t = zeros(s.d,1);
p.robInd_R = eye(s.d);
if ( s.gt_bool==1 ) % translate estimate if ground truth is available
    p.robInd_t = rob(robInd).t(:,1);
    p.robInd_R = rob(robInd).R(:,:,1);
else
    p.robInd_t = rob(robInd).t_hat_d(:,1);
    ind_robInd = (1:s.d)+(robInd-1)*s.d;
    p.robInd_R = rob(robInd).R_hat_d(:,:,1);
end

% build a spanning tree through via bredth first search
% to visualize distributed estimate results
rootInd = robInd;
if ( ~exist('s.G','var') )
    s = buildSpanningTree(s,rootInd);
end
[p.rob.t_plot,p.rob.R_plot] = spanningTreeTransformations(s,rob,rootInd,tInd,0);

% robot's estimate to plot
i = robInd;
ind_i = (1:s.d)+(i-1)*s.d;
% transform estimation about frame i in world frame for comparison
t_trans = p.robInd_t - p.rob.t_plot(:,i);
R_trans = p.rob.R_plot(:,ind_i)'*p.robInd_R;
[p.t_cur, p.R_cur] = ...
    rigidBodyTrans(p.rob.t_plot,p.rob.R_plot,t_trans,R_trans,i);
%    p.t_cur = p.rob.t_plot;
%    p.R_cur = p.rob.R_plot;
% remove agents with no estimate from ith point of view
zeroInd = find(p.rob.t_plot(1,:)==0); % indices to remove from lack of relative measurement
zeroInd(zeroInd==i) = [];
% plot communication network
for e = 1:s.m
    i = s.M(e,1);
    j = s.M(e,2);
    x_pts = [p.t_cur(1,i); p.t_cur(1,j)];
    y_pts = [p.t_cur(2,i); p.t_cur(2,j)];
    if ( s.d==2 )
        plot(x_pts,y_pts, ...
            'Color',p.commEstColor,'LineStyle','-','LineWidth',1);
    else
        z_pts = [p.t_cur(3,i); p.t_cur(3,j)];
        plot3(x_pts,y_pts,z_pts, ...
            'Color',p.commEstColor,'LineStyle','-');
    end
end
% plot coordinate systems
if ( s.n<p.maxPlotN )
    for ii = 1:s.n
        if ( any(p.t_cur(:,ii)~=0) && ~any(ii==zeroInd) )
            ind_ii = (1:s.d)+(ii-1)*s.d;
            % for coordinate system labels
            % plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
            %              ['$\hat{',num2str(ii),'}$'], 50, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
            % no coordinate system labels
            plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
                '', 0, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
        end
    end
end
hold off;
box on;
rotate3d on;
axis equal;
view(p.plotViewValues);
title(titleString,'Interpreter','LaTex');
% define plot volume
p.plotVol = zeros(1,6);
min_bound = 1E6*ones(3,1);
max_bound = -1E6*ones(3,1);
for i = 1:s.n
    ind = find(p.t_cur(:,i)<min_bound);
    for j = ind'
        p.plotVol((j-1)*2+1) = p.t_cur(j,i);
        min_bound(j) = p.t_cur(j,i);
    end
    ind = find(p.t_cur(:,i)>max_bound);
    for j = ind'
        p.plotVol((j-1)*2+2) = p.t_cur(j,i);
        max_bound(j) = p.t_cur(j,i);
    end
end
offset = 0.1;
p.plotVol = p.plotVol + [-offset,offset,-offset,offset,-offset,offset];
axis(p.plotVol);
% axis off; % to hide axes and background



%% Plot Distributed Result for Multi-UAV Experiment

if strcmp(init.pose, 'experiment')
    
    % Initialize Figure
    figure(4);
    clf(4);
    titleString = 'GeoD Result with 3D Reconstruction';
    hold on;
    im_indices = [194, 4]; % define images to be displayed in figure window
    robInd = 165; % plot estimate from robInd's perspective (165)
    p.robInd_t = zeros(s.d,1);
    p.robInd_R = eye(s.d);
    p.robInd_t = rob(robInd).t(:,1);
    p.robInd_R = rob(robInd).R(:,:,1);

    % toggle diffenent experimental results
    plot_result = 'geod'; % geod, sesync, or ground_truth
    tInd_ = length(mle_total); % convergence for distributed
    for i = 1:s.n
        ind_i = (1:s.d)+(i-1)*s.d;
        if strcmp(plot_result, 'geod')
            p.rob.t_plot(:,i) = rob(i).t_hat_d(:,tInd_);
            p.rob.R_plot(:,ind_i) = rob(i).R_hat_d(:,:,tInd_);
        elseif strcmp(plot_result, 'sesync')
            p.rob.t_plot(:,i) = rob(i).t_hat(:,1);
            p.rob.R_plot(:,ind_i) = rob(i).R_hat(:,:,1);
        elseif strcmp(plot_result, 'ground_truth')
            p.rob.t_plot(:,i) = rob(i).t(:,1);
            p.rob.R_plot(:,ind_i) = rob(i).R(:,:,1);
        else
            error('Error in plotResults.m: incorrect type of result provided');
        end
    end
    
    % robot's estimate to plot
    i = robInd;
    ind_i = (1:s.d)+(i-1)*s.d;
    % transform estimation about frame i in world frame for comparison
    t_trans = p.robInd_t - p.rob.t_plot(:,i);
    R_trans = p.rob.R_plot(:,ind_i)'*p.robInd_R;
    [p.t_cur, p.R_cur] = ...
        rigidBodyTrans(p.rob.t_plot,p.rob.R_plot,t_trans,R_trans,i);
    % remove agents with no estimate from ith point of view
    zeroInd = find(p.rob.t_plot(1,:)==0); % indices to remove from lack of relative measurement
    zeroInd(zeroInd==i) = [];
    % plot coordinate systems
    if ( s.n<p.maxPlotN )
        for ii = 1:s.n
            if ( any(p.t_cur(:,ii)~=0) && ~any(ii==zeroInd) )
                ind_ii = (1:s.d)+(ii-1)*s.d;
                if any(im_indices==ii)
                    plotCameraSimple([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii)*R_rc)],...
                        [0 640 0 480],-1, p.cVec(p.cInd(ii),:),1.0*p.arrowLength,1,1);
                else
                    % plot camera frames
                    plotCameraSimple([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii)*R_rc)],...
                        [0 640 0 480],-1, p.cVec(p.cInd(ii),:),1.0*p.arrowLength,1,0);
                    % % for coordinate system labels
                    % plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
                    %              '', 0, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
                    % % % no coordinate system labels
                    % plotCoordSys([p.t_cur(:,ii);rot2euler(p.R_cur(:,ind_ii))], ...
                    %              ['$\hat{',num2str(ii),'}$'], 50, p.cVec(p.cInd(ii),:), p.estAlpha, 1.0*p.arrowLength, 2);
                end
            end
        end
    end
    % plot communication network
    for e = 1:s.m
        i = s.M(e,1);
        j = s.M(e,2);
        x_pts = [p.t_cur(1,i); p.t_cur(1,j)];
        y_pts = [p.t_cur(2,i); p.t_cur(2,j)];
        z_pts = [p.t_cur(3,i); p.t_cur(3,j)];
        plot3(x_pts,y_pts,z_pts, ...
            'Color',p.commColor,'LineStyle','-','LineWidth',1);
    end
    
    % display 3d model
    if ( strcmp(init.pose,'experiment') )
        if ~exist('object','var')
            object = pcread([experiment_data_path,'filtered_dense_cloud.ply']);
        end
        % scale point cloud
        ptCldScale = 0.085;
        A = [ptCldScale 0 0 0; ...
            0 ptCldScale 0 0; ...
            0 0 ptCldScale 0; ...
            0 0 0 1];
        tform = affine3d(A);
        scaledObject = pctransform(object,tform);
        % transform point cloud
        A = [euler2rot([-pi/2+0.04; -0.24; -0.06]), [0; 0; 0];...
            zeros(1,3),1];
        tform = affine3d(A);
        transformedObject = pctransform(scaledObject,tform);
        D = zeros(size(transformedObject.Location));
        D(:,1) = 0.2;
        D(:,2) = -2.2;
        D(:,3) = -min(transformedObject.ZLimits);
        transformedObject = pctransform(transformedObject,D);
        pcshow(transformedObject);
    end
    
    % axes and view
    hold off;
    rotate3d on;
    box off;
    title(titleString,'Interpreter','LaTex');
    axis equal;
    view([135, 25]); % iso
    % view([145, 90]); % top
    axis off; % to hide axes and background
    
end


%% Distributed Optimization Plots

robInd = 1;
i = robInd;
ind_i = (1:s.d)+(i-1)*s.d;
totalSteps_long = s.steps;
totalSteps = s.steps;
% Partial Simulation
if (tInd~=s.steps)
    totalSteps = tInd;
end

% objective value plot
figure(10); clf(10);
% translation term
subplot(4,1,1); hold on; box on;
title('Pose Graph Objective Value','Interpreter','LaTex');
if sesync_comparison
    plot(1:totalSteps,repmat(s.mle_error(1,1),1,totalSteps),'r','LineWidth',2);
end
plot(1:totalSteps,s.mle_error_d(1,1:totalSteps),'r-.','LineWidth',2);
hold off;
if sesync_comparison
    legend('SE-Sync','GeoD');
end
ylabel('Translation Value','Interpreter','LaTex');
axis([1, max(totalSteps, 2), 0, max([s.mle_error_d(1,1:totalSteps),1])]);
% chordal rotation term
subplot(4,1,2); hold on; box on;
if sesync_comparison
    plot(1:totalSteps,repmat(s.mle_error(2,1),1,totalSteps),'g','LineWidth',2);
end
plot(1:totalSteps,s.mle_error_d(2,1:totalSteps),'g-.','LineWidth',2);
hold off;
if sesync_comparison
    legend('SE-Sync','GeoD');
end
ylabel('Chordal Value','Interpreter','LaTex');
axis([1, max(totalSteps, 2), 0, max([s.mle_error_d(2,1:totalSteps),1])]);
% geodesic rotation term
subplot(4,1,3); hold on; box on;
if sesync_comparison
    plot(1:totalSteps,repmat(s.mle_error(3,1),1,totalSteps),'b','LineWidth',2);
end
plot(1:totalSteps,s.mle_error_d(3,1:totalSteps),'b-.','LineWidth',2);
hold off;
if sesync_comparison
    legend('SE-Sync','GeoD');
end
ylabel('Geodesic Value','Interpreter','LaTex');
axis([1, max(totalSteps, 2), 0, max([s.mle_error_d(3,1:totalSteps),1])]);
hold off;
xlabel('Iteration','Interpreter','LaTex');
% total objective value plot
subplot(4,1,4); hold on; box on;
% % chordal objective
% plot(1:totalSteps,repmat(s.mle_error(1,1)+s.mle_error(2,1),1,totalSteps),'b','LineWidth',2);
% plot(1:totalSteps,s.mle_error_d(1,1:totalSteps)+s.mle_error_d(2,1:totalSteps),'b-.','LineWidth',2);
% geodesic objective
if sesync_comparison
    plot(1:totalSteps,repmat(s.mle_error(1,1)+s.mle_error(3,1),1,totalSteps),'b','LineWidth',2);
end
plot(1:totalSteps,s.mle_error_d(1,1:totalSteps)+s.mle_error_d(3,1:totalSteps),'b-.','LineWidth',2);
axis([1, max(totalSteps, 2), 0, max([s.mle_error_d(1,1:totalSteps)+s.mle_error_d(3,1:totalSteps),1])]);
hold off;
xlabel('Iteration','Interpreter','LaTex');
ylabel('Total Geodesic Value','Interpreter','LaTex');
if sesync_comparison
    legend('SE-Sync','GeoD');
end

