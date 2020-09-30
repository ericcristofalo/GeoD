%--------------------------------------------------------------------------
%
% File Name:      randRot.m
% Date Created:   2018/02/04
% Date Modified:  2018/02/04
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Generate a random rotation matrix in SO(d)
%                 Method: Arvo, James. "Fast random rotation matrices." 
%                         Graphics Gems III (IBM Version). 1992. 117-120.
%
% Inputs:         d: d in {2,3} 
%                 n: number of desired rotation matrices (optional)
%
% Outputs:        R: R in SO(d)
%                    if n = null or 1: R has size [d,d]
%                    if n >1: R has size [d,d*n]
%
% Example:        R = randRot(3)
%
%--------------------------------------------------------------------------

function R = randRot(varargin)

% Handle Input Arguments
if (nargin==0)
   d = 3;
   n = 1;
elseif (nargin==1)
   d = varargin{1};
   n = 1;
elseif (nargin==2)
   d = varargin{1};
   n = varargin{2};
else
   error('Error in randRot.m: too many input arguments');
end

% Error Handling
if (d~=2 && d~=3)
   error('Error in randRot.m: dimension must be either 2 or 3');
elseif (n<1)
   error('Error in randRot.m: number of rotation matrices must be >=1');
end

R = zeros(d,n*d);
x = rand(d,n);
for i = 1:n
   ind_i = (1:d)+(i-1)*d;
   if (d==2)
      R(:,ind_i) = rotMat(2*pi*rand(1,1)-pi);
   else
      v = [cos(2*pi*x(2,i))*sqrt(x(3,i));
      sin(2*pi*x(2,i))*sqrt(x(3,i));
      sqrt(1-x(3,i))];
      H = eye(d)-2*(v*v');
      R(:,ind_i) = -H*rotMat(2*pi*x(1,i),'z');
   end
end

end