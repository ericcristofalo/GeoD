%--------------------------------------------------------------------------
%
% File Name:      expMap.m
% Date Created:   2017/09/21
% Date Modified:  2018/06/25
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Exponential Map for Rotation Rates
%                 Tron et al. 
%                 Distributed Pose Averaging in Camera Networks via Consensus on SE(3)
%                 Made some changes for robustness
%
% Inputs:         omega: vector of rotation rates, either 1 or 3
%                        dimensional (d=2 or d=3)
%
% Outputs:        R: exp(omega)
%
%--------------------------------------------------------------------------

function R = expMap(omega)

% Check InputSize
if ( size(omega,2)~=1 )
  omega = vee(omega);
end
d_ = length(omega);

% Compute Skew Symmetric Matrix
if ( d_==1 )
  d = 2;
  Omega = [0 -omega; omega 0];
elseif ( d_==3 )
  d = 3;
  Omega = hat(omega);
else
  error('Error in expMap.m: wrong dimension of rotation rate vector');
end

% Compute Theta
theta = norm(omega);
vector = omega;
if theta > 2*pi
   theta2 = theta;
   theta = mod(theta,2*pi);
   vector = vector*theta/theta2;
end
if theta > pi
   vector = vector*(1 - 2*pi/theta);
end
theta = norm(vector);

% Compute Rotation Matrix
normOmega = norm(omega);
if ( normOmega==0 )
  R = eye(d);
else
  R = eye(d) + ...
      Omega/normOmega*sin(theta) + ...
      Omega*Omega/(normOmega^(2.0))*(1.0-cos(theta));
end

end