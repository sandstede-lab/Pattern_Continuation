% clear all MATLAB variables, command windows, and figure
function [u,v] = Schnakenberg(a,b, viz)

if nargin < 1
    a = 2.5;
end
if nargin < 2
    b = 3;
end
if nargin < 3
    viz = 0;
end
D1 = 0.005;
D2 = 1;



% spatial discretization
xright= 4; nx=80;%nx = 50;
% yright=128; ny=nx;
% xright = 10; nx = 100;
dx=xright/nx;




% set the parameters (time discretization)
dt=0.02*dx^2;
dt = 0.001*2;
dt = 0.001*0.5;


maxit=round(1.5*500/dt);
% maxit = round(100/dt);
nn=maxit;

% for pit=1:50    % the number of images
% set the initial condition

u = b + zeros(nx);
v = (b+a)/2/b/b + 0.1*rand(nx,nx);
% r = x + 0.001*rand(nx,nx);

for t = 1:nn
    up = u; vp = v; %rp = r;
    
    f = - u + u.^2.*v + (b-a)/2;
    g = -u.^2.*v + (b+a)/2;
%     h = 1/2/epsilon*(x-r);
    
    u = dt*(D1 * laplacian(up, dx) +  f) + up;   % Main Iteration Loops for u and v
    
    v = dt*(D2 * laplacian(vp, dx) + g) + vp;
    
%     r = dt*(Dr * laplacian(rp, dx) + h) + rp;
    time_elapsed = t*dt;  % time progress to be displayed on the screen as the program runs
%     t_max
%     nt
%     t

% Display of results

% mesh(x,y,u) % creates a mesh plot of the results

if viz && mod(t, round(1/dt)) == 0
    figure(1);
    % [X,Y] = meshgrid(x,y);
    imagesc(u)
    colorbar()
    title(num2str(dt*t))
    pause(0.1)
end

end