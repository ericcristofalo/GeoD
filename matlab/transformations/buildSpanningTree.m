%--------------------------------------------------------------------------
%
% File Name:      buildSpanningTree.m
% Date Created:   2018/07/27
% Date Modified:  2018/07/27
%
% Author:         Eric Cristofalo
% Contact:        eric.cristofalo@gmail.com
%
% Description:    Function for building a spanning tree through a connected graph
%
% Inputs:
%
% Outputs:
%
%--------------------------------------------------------------------------

function s = buildSpanningTree(s,rootInd)

   s.G = cell(s.n,2); % cell of parents and children nodes {{parents},{childen}}
   robList = 1:s.n; % list of node indices starting with the rootInd
   nodeList = [];
   
   while ( ~isempty(robList) )
      
      if ( any(robList==rootInd) ) % begin with root of tree
         i = rootInd;
         robList(robList==i) = [];
      else % choose next node from list
         % if nodeList is empty before robList, graph is disjoint
         if ( isempty(nodeList) )
            disp('oh no!')
         end
         i = nodeList(1,1); % choose first node added
         nodeList = nodeList(2:end); % remove from nodeList
      end
      % Find Neighbor Indices
      ind_m = find(s.M(:,1)==i);
      ind_m = s.M(ind_m,:);
      for j = ind_m(:,2)'
         % If Child Not in Tree, Add to Tree
         if ( any(robList==j) )
            s.G{i,2} = [s.G{i,2},j]; % add child index
            s.G{j,1} = [s.G{j,1},i]; % add parent index
            robList(robList==j) = []; % remove from robList
            nodeList = [nodeList,j]; % add to nodeList
         end
      end
      
   end
   
end