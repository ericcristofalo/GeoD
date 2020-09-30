%--------------------------------------------------------------------------
%
% File Name:      plotCameraSimple.m
% Date Created:   2014/11/15
% Date Modified:  2014/11/15
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Plots cartesian coordinate system given position w.r.t a
%                 global coordinate system.
%
% Inputs:         pose: vector with the entires [x; y; z; phi; theta; psi]
%                       Pose is w.r.t. the world frame, angles measured in 
%                       radians.
%                 plotSize: virtual size of image [0 640 0 480]
%                 lineColor: vector denoting the color of the lines from 
%                            camera center to virtual image frame.
%                            set to -1 for no lines drawn
%                 imageColor: vector denoting the color of virtual image
%                 arrowLength: desired arrow length.
%                 lineThick: desired line thickness.
%                 filledIn: set to 1 if image should be filled
%
% Outputs:        plot
%
% Example:        plotCamera([1; 1; -1; 0.1; 0.1; 0.7],[0 640 0 480],[0,0,1],[1,1,1],0.5,2,0);
%
%--------------------------------------------------------------------------

function plotCameraSimple(pose, plotSize, lineColor, imageColor, arrowLength, lineThick, filledIn)

% Plot Camera
R_wc = euler2rot(pose(4:6));
scale = arrowLength;
ratio = (plotSize(4)-plotSize(3))/...
         (plotSize(2)-plotSize(1)); % height/width
depth = arrowLength;
if ( lineColor(1)==-1 )
    depth = 0;
end
corners = [scale/2,scale/2*ratio,depth;...
           -scale/2,scale/2*ratio,depth;...
           -scale/2,-scale/2*ratio,depth;...
           scale/2,-scale/2*ratio,depth];
for i = 1:4
   dirVec = corners(i,:)'; % current axis direction vector
   dirVec = R_wc*dirVec;
   if ( lineColor(1)~=-1 && arrowLength~=0.0 )
       quiver3(pose(1),pose(2),pose(3),...
               dirVec(1),dirVec(2),dirVec(3),...
               1,'Color',lineColor);
   end
   corners(i,:) = pose(1:3)'+dirVec';
end
corners = [corners; corners(1,:)];
if filledIn
   fill3(corners(:,1),corners(:,2),corners(:,3),imageColor);
else
   plot3(corners(:,1),corners(:,2),corners(:,3),'Color',imageColor,'LineWidth',lineThick);
end


end