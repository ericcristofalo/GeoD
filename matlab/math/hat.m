%--------------------------------------------------------------------------
%
% File Name:      hat.m
% Date Created:   2014/08/06
% Date Modified:  2019/03/19
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Hat Operator
%                 Creates a skew-symmetric matrix given a vector
%
% Inputs:         v: 1D or 3D column vector
%
% Outputs:        S: [2x2] or [3x3] skew-symmetric matrix
%
% Example:        v = [1;2;3]
%                 S = skewSymMat(v)
%                 S =
%                     0    -3     2
%                     3     0    -1
%                    -2     1     0
%
%--------------------------------------------------------------------------

function S = hat(v)
   if size(v,1)==1
      S = [ 0  -v;
            v   0];
   elseif size(v,1)==3
      S = [    0  -v(3)   v(2) ;
            v(3)      0  -v(1) ;
           -v(2)   v(1)      0 ];
   else
      error('Incorrect size of matrix input to hat.m');
   end
end