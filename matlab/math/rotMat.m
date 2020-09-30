%--------------------------------------------------------------------------
%
% File Name:      rotMat.m
% Date Created:   2014/04/07
% Date Modified:  2017/09/19
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Rotation matrix calculaton function.
%
% Inputs:         angle: angle of rotation (in radians)
%                 axis: axis about which the rotation takes place (string)
%
% Outputs:        R: the resulting rotation matrix
%
% Example:        R = rotMat(pi/3);
%                 R = rotMat(psi,'z')*rotMat(the,'y')*rotMat(phi,'x');
%
%--------------------------------------------------------------------------

function  R = rotMat(varargin)

if nargin==1
   % Planar Rotation Matrix
   angle = varargin{1};
   R = [cos(angle)  -sin(angle);
        sin(angle)   cos(angle)];
elseif nargin==2
   % 3D Rotation Matrix
   if ischar(varargin{1})
      axis = varargin{1};
      angle = varargin{2};
   else
      angle = varargin{1};
      axis = varargin{2};
   end
   switch axis
      case 'x'
         R = [1           0           0          ;
              0           cos(angle)  -sin(angle);
              0           sin(angle)  cos(angle) ];
      case 'y'
         R = [cos(angle)  0           sin(angle) ;
              0           1           0          ;
              -sin(angle) 0           cos(angle) ];
      case 'z'
         R = [cos(angle)  -sin(angle) 0          ;
              sin(angle)  cos(angle)  0          ;
              0           0           1          ];
   end
else
  error('Error in rotMat.m: wrong number of inputs');
end

end