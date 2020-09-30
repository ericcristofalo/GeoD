function rob = initRobDataStructure(s)

rob(s.n) = struct();
% robot world frame initialization
for i = 1:s.n
    % ground truth only for certain datasets
    if ( s.gt_bool )
        rob(i).t = zeros(s.d,1);     % world position
        rob(i).R = zeros(s.d,s.d,1); % row vector or rotation matrices
    end
    % measurements
    rob(i).y_t = zeros(s.d,s.n,1);     % measurements
    rob(i).y_R = zeros(s.d,s.d*s.n,1); % measurements
    % centralized estimates
    rob(i).t_hat = zeros(s.d,1);     % translation
    rob(i).R_hat = zeros(s.d,s.d,1); % rotation
end

end

