%--------------------------------------------------------------------------
%
% File Name:      euler2rot.m
% Date Created:   2014/06/27
% Date Modified:  2017/09/19
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Rotation matrix computed from 3 Euler angles.
%
% Inputs:         euler: [3x1] vector of Euler angles or angle about z-axis
%                        in 2D system
%
% Outputs:        R: [3x3] rotation matrix or [2x2] rotation matrix
%
% Example:        R = euler2rot(euler)
%
%--------------------------------------------------------------------------

function R = euler2rot(euler)

if ( size(euler,1)==1 )
  R = rotMat(euler);
elseif ( size(euler,1)==3 )
  phi = euler(1);
  the = euler(2);
  psi = euler(3);
  R = rotMat(psi,'z')*rotMat(the,'y')*rotMat(phi,'x');
else
  error('Error in euler2rot.m: wrong number of Euler angles used.');
end

end