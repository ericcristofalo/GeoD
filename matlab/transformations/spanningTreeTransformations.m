%--------------------------------------------------------------------------
%
% File Name:      spanningTreeTransformations.m
% Date Created:   2018/07/27
% Date Modified:  2020/01/22
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Searches a spanning tree for each node back to a 
%                 desired root node
%
% Inputs:
%
% Outputs:
%
%--------------------------------------------------------------------------

function [t_tree, R_tree] = spanningTreeTransformations(s,rob,rootInd,tInd,useMeasurements)
   
   % Initialization
   t_tree = zeros(s.d,s.n); % array of all translation vectors
%    R_tree = zeros(s.d,s.d*s.n); % array of all stacked rotation matrices
   R_tree = repmat(eye(s.d), 1, s.n);
   
   % Define Search Order Based on rootInd Node
   robList = [1:rootInd-1,rootInd+1:s.n]; % list of node indices starting with the rootInd
   robList = [rootInd,robList];
   
   % For Each Robot
   for i = robList
      ind_i = (1:s.d)+(i-1)*s.d;
      % Search Parents Back to rootInd Root Node
      if ( i==rootInd ) % begin with root of tree
         t_tree(:,1) = zeros(s.d,1); % initialize root to identity transform since tree is relative!
         R_tree(:,1:s.d) = eye(s.d);
      else % choose random node from list
         robTemp = i;
         t_temp = zeros(s.d,1);
         R_temp = eye(s.d);
         while ( robTemp~=rootInd ) % recursively compute transforms back to root of tree
            parNode = s.G{robTemp,1}; % find parent of node i
            ind_par = (1:s.d)+(parNode-1)*s.d;
            ind_rob = (1:s.d)+(robTemp-1)*s.d;
%             % Centralized Testing Only:
%             R_cur = rob(parNode).R_hat(:,ind_temp,1); % R_par_i
%             t_temp = rob(parNode).t_hat(:,robTemp,1) + ...
%                      R_cur*t_temp; % estimate of t_par_i
            % Set Relative Transformation From Estimates (can be in any frame now)
            if ( useMeasurements )
               R_cur = ...
                  rob(parNode).y_R(:,ind_rob,tInd); % R_par_rob = R_w_par'*R_w_rob
               t_temp = ...
                  rob(parNode).y_t(:,robTemp,tInd) + ...
                  + R_cur*t_temp; % estimate of t_par_rob = t_par_rob + R_par_rob*t_rob_child
            else
               R_cur = ...
                  rob(parNode).R_hat_d(:,:,tInd)'*...
                  rob(robTemp).R_hat_d(:,:,tInd); % R_par_rob = R_w_par'*R_w_rob
               t_temp = ...
                  rob(parNode).R_hat_d(:,:,tInd)'*...
                  (rob(robTemp).t_hat_d(:,tInd) - ...
                  rob(parNode).t_hat_d(:,tInd)) ...
                  + R_cur*t_temp; % estimate of t_par_rob = t_par_rob + R_par_rob*t_rob_child
            end
            R_temp = R_cur*R_temp;
            robTemp = parNode;
         end
         t_tree(:,i) = t_temp;
         R_tree(:,ind_i) = R_temp;
      end
   end
   
end