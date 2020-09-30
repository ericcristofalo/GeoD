%--------------------------------------------------------------------------
%
% File Name:      plotCoordSys.m
% Date Created:   2014/11/15
% Date Modified:  2017/08/04
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Plots cartesian coordinate system given position w.r.t a
%                 global coordinate system.
%                 All inputs in meters. 
%
% Inputs:         pose: vector with the entires [x; y; z; phi; theta; psi]
%                       Pose is w.r.t. the world frame, angles are measured  
%                       in radians.
%                 name: robot name for text plot if desired
%                 diameter: diameter of centroid sphere
%                 color: vector denoting color of centroid sphere [0,0,0]
%                 alpha: scale for orthogonal line colors
%                 arrowLength: desired arrow length.
%                 lineThick: desired line thickness.
%
% Outputs:        the plot
%
% Example:        plotCoordSys([1; 1; -1; 0.1; 0.1; 0.7],1,100,[1,0,0],1,0.5,1);
%
%--------------------------------------------------------------------------

function [h] = plotCoordSys(pose, name, diameter , color, alpha, arrowLength, lineThick)

% Size of Problem
if ( size(pose,1)==3 )
  d = 2;
  r_ind = 3;
elseif ( size(pose,1)==6 )
  d = 3;
  r_ind = 4:6;
else
  error('Error in plotCoordSys.m: wrong size of pose vector. Must be either length 3 or length 6.');
end

% Plot Coordinate System Orientation
R_wr = euler2rot(pose(r_ind));
vecColor = [1,0,0,alpha; 0,1,0,alpha; 0,0,1,alpha];
for i = 1:d % for each direction (x,y,z)
  dirVec = zeros(d,1); % current axis direction vector
  dirVec(i) = arrowLength;
  dirVec = R_wr*dirVec;
  if ( d==2 )
%     quiver(pose(1),pose(2),...
%            dirVec(1),dirVec(2),...
%            arrowLength,'color',vecColor(i,:),...
%            'LineWidth',lineThick);
    plot([pose(1),pose(1)+dirVec(1)],...
         [pose(2),pose(2)+dirVec(2)],...
         'color',vecColor(i,:),...
         'LineWidth',lineThick);
  elseif ( d==3 )
%     quiver3(pose(1),pose(2),pose(3),...
%             dirVec(1),dirVec(2),dirVec(3),...
%             arrowLength,'color',vecColor(i,:),...
%             'LineWidth',lineThick);
    plot3([pose(1),pose(1)+dirVec(1)],...
          [pose(2),pose(2)+dirVec(2)],...
          [pose(3),pose(3)+dirVec(3)],...
          'color',vecColor(i,:),...
          'LineWidth',lineThick);
  end
end

% Plot Coordinate System Position
if ( d==2 )
   if ( diameter~=0 )
      scatter(pose(1,1),pose(2,1),diameter,color,'Fill');
   end
elseif ( d==3 )
   if ( diameter~=0 )
      scatter3(pose(1,1),pose(2,1),pose(3,1),diameter,color,'Fill');
   end
end

% Plot Coordinate System ID Text
offset = arrowLength/2.0;
% tempPose = R_wr*offset*ones(3,1);
tempPose = offset*[1;0;1];
if ( d==2 )
  h = text(pose(1)+tempPose(1),...
           pose(2)+tempPose(2),...
  name,'HorizontalAlignment','left','FontSize',12,'Interpreter','LaTex');
elseif (d==3)
  h = text(pose(1)+tempPose(1),...
           pose(2)+tempPose(2),...
           pose(3)+tempPose(3),...
  name,'HorizontalAlignment','left','FontSize',12,'Interpreter','LaTex');
end

end

