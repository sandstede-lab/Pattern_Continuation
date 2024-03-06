function [ lap ] = laplacian( z, hx )
% This function calculates the laplacian of a given variable Z stored as a
% matrix z(i,j)and where z(i,j) represents value of Z at (i,j)th point on
% the spatial grid. The boundary conditions are periodic.

n = size(z, 2);

% lap = -4*z + z(:,[2:n,1]) + z(:,[n,1:n-1]) + z([2:n,1],:) + z([n,1:n-1],:);

lap = -4*z + (circshift(z,1)+circshift(z,-1)+circshift(z,-1,2)+circshift(z,1,2));

lap = lap/hx/hx;                    
   

            
        
        



end

