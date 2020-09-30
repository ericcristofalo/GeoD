%--------------------------------------------------------------------------
%
% File Name:      rigidBodyTrans
% Date Created:   2017/08/15
% Date Modified:  2018/03/04
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Applies rigid body transformation (translation,rotation) 
%                 to a set of 3D poses given in the world frame. The set of
%                 poses denotes a rigid body in the world frame. 
%
% Inputs:         t: [d x n] matrix of doubles is the list of positions
%                   meansure in the world frame, where d is the dimension
%                   of the pose and n is the number of poses. 
%                 R: [d x d*n] matrix of doubles is the row vector of
%                   concatonated individual rotation matrices in SO(d) 
%                   with respect to the world frame. 
%                 translation: [d x 1] vector of doubles denotes the
%                   rigid body translation measured in world frame.
%                 rotation: [d x d] matrix of doubles denotes the
%                   rigid body rotation applied to the rigid body about 
%                   the specified point with respect to the world frame.
%                 point: [d x 1] vector for the point in world frame 
%                 	(before translation applied) about which to make the 
%                   rotation.
%                   Setting point = 0 indicates: point = COM of rigid body.
%                   Setting point = i indicates: point = ith pose in set. 
%
% Outputs:        t_new: [d x n] matrix of doubles is the new positions 
%                   in the world frame.
%                 R_new: [d x d*n] matrix of doubles is the individual 
%                   rotation matrices with respect to the world frame.
%
% Example:        see rigidBodyTrans_Testing.m script. 
%
%--------------------------------------------------------------------------

function  [t_new, R_new] = rigidBodyTrans(t,R,translation,rotation,point)

% Dimensions
d = size(t,1);
n = size(t,2);

% Transformation Frame
% The frame about which the rotation is made
t_origin = zeros(d,1);
R_w_origin = eye(d);
if ( size(point,1)==1 ) % use the point from the set of poses
  if point==0 % use the COM of the set of poses
    for i = 1:n
      t_origin = t_origin + t(:,i);
    end
    t_origin = t_origin./n;
  else % use the ith pose in the set of poses
    if ( point<=size(t,2) )
      t_origin = t(:,point);
      ind_frame = (1:d)+(point-1)*d;
      R_w_origin = R(:,ind_frame);
    else
      error('Error in rigidBodyTrans.m: point exceeds dimensions of input t matrix.')
    end
  end
elseif ( size(point,1)==3 ) % use specifed point in world frame
  t_origin = point;
else
  error('Error in rigidBodyTrans.m: wrong size of point.');
end

% Perform Rigid Body Transformation
R_new = zeros(d,d*n);
t_new = zeros(d,n);
for i = 1:n
  ind_i = (1:d)+(i-1)*d;
  R_new(:,ind_i) = R_w_origin*rotation*R_w_origin'*R(:,ind_i);
  t_new(:,i) = t_origin + translation + ...
               R_w_origin*rotation*R_w_origin'*(t(:,i)-t_origin);
end

% % Perform Rigid Body Transformation
% R_new = zeros(d,d*n);
% t_new = zeros(d,n);
% R_new(:,ind_frame) = R_w_origin*rotation; % set origin equal to (t_origin, R_w_origin)
% t_new(:,frame) = t_origin + translation;
% for i = [1:frame-1,frame+1:n]
%   ind_i = [i*d-2, i*d-1, i*d-0];
%   R_new(:,ind_i) = R_new(:,ind_frame)*R_w_origin'*R(:,ind_i);
%   t_new(:,i) = t_new(:,frame) + R_new(:,ind_frame)*R_w_origin'*(t(:,i)-t_origin);
% end

end

