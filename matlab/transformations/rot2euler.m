%--------------------------------------------------------------------------
%
% File Name:      rot2euler.m
% Date Created:   2014/07/24
% Date Modified:  2018/03/23
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Calculates the Euler angles given a rotation matrix.
%
% Inputs:         R: rotation matrix [dxd] or [dx(nd)] where n is the
%                 number of rotation matrices in a rotation matrix. 
%
% Outputs:        euler: vector of Euler angles [3xn] or [n]
%
% Example:        euler = rot2euler(R)
%
%--------------------------------------------------------------------------

function euler = rot2euler(R)

d1 = size(R,1);
d2 = size(R,2);

if ( d1==2 )
   if ( d2==2 )
      euler = acos(R(1,1));
   elseif ( mod(d2,2)==0 )
      euler = zeros(1,d2);
      for i = 1:(d2/2)
         euler(1,i) = acos(R(1,1+(i-1)*2));
      end
   else
      error('Error in rot2euler.m: wrong dimensions of input rotation matrix');
   end
elseif ( size(R,1)==3 )
   if ( d2==3 )
      euler = zeros(3,1);
      euler(1,1) = atan2(R(3,2),R(3,3));	% phi
      euler(2,1) = asin(-R(3,1));      	% theta
      euler(3,1) = atan2(R(2,1),R(1,1)); 	% psi
   elseif ( mod(d2,3)==0 )
      euler = zeros(3,d2/3);
      for i = 1:(d2/3)
         R_temp = R(:,(1:3)+(i-1)*3);
         euler(1,i) = atan2(R_temp(3,2),R_temp(3,3)); % phi
         euler(2,i) = asin(-R_temp(3,1));             % theta
         euler(3,i) = atan2(R_temp(2,1),R_temp(1,1)); % psi
      end
   else
      error('Error in rot2euler.m: wrong dimensions of input rotation matrix');
   end
   
else
   error('Error in rot2euler.m: wrong dimensions of input rotation matrix');
end

end