%--------------------------------------------------------------------------
%
% File Name:      vee.m
% Date Created:   2014/08/06
% Date Modified:  2019/03/19
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Vee Operator
%                 Decomposes a skew-symmetric matrix into a vector
%
% Inputs:         S: [2x2] or [3x3] skew-symmetric matrix
%
% Outputs:        v: 1D or 3D column vector
%
% Example:        S = [0,-3,2; 3,0,-1; -2,1,0]
%                 v = skewSymMatInv(S)
%                 v =
%                     1
%                     2
%                     3
%
%--------------------------------------------------------------------------

function v = vee(S)
   if size(S,1)==2
      v = S(2,1);
   elseif size(S,1)==3
      v = [-S(2,3);S(1,3);-S(1,2)];
   else
      error('Incorrect size of vector input to skewSymMatInv.m');
   end
end