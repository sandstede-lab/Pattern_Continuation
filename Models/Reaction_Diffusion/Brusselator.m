function [u,v,x,y] = Brusselator(a,b,viz,Nx,t_max,Uinit, Vinit)

% Reaction Diffusion Equations:Brusselator Equations
 
% Initialization



Du = 4;    
Dv = 32;
if nargin < 3
    viz = 0;
end
if nargin < 4
    Nx = 50;    % spatial grid cell length
end
ht = 0.005;    % time step
% ht = 0.01;
if nargin < 5
    t_max = 2*100;
end
% t_max = 2*t_max;
x_max = 50;
x = linspace(0, x_max, Nx);
hx = x(2) - x(1);
ht = ht*hx^2;
y = x;
n = size(x, 2);  % number of time iterations
nt = t_max/ht; % number of space points along one direction
u = zeros(n,n); % initialization of quantity u and v
v = u;

brusselator_f = @(p,q,a,b) a - (b+1)*p + q.*p.^2;
brusselator_g = @(p,q,a,b) b*p - q.*p.^2;

if nargin < 6 || size(Uinit,1) == 0
    u = a + 0.05*randn(n,n);
    v = v + b/a;
else
    u = Uinit;
    v = Vinit;
end

% Main Loop

for t = 1:nt
    p = u;
    q = v;
    
    u = ht*(Du*laplacian(p, hx) + brusselator_f(p, q, a, b)) + p;   % Main Iteration Loops for u and v
    v = ht*(Dv*laplacian(q, hx) + brusselator_g(p, q, a, b)) + q;   
    
    time_elapsed = t*ht;  % time progress to be displayed on the screen as the program runs
%     t_max
%     nt
%     t

% Display of results

% mesh(x,y,u) % creates a mesh plot of the results

if viz && mod(t, round(1/ht)) == 0
    figure(1);
    % [X,Y] = meshgrid(x,y);
    imagesc(u)
    colorbar()
    title(num2str(ht*t))
    pause(0.1)
end

            
            
            
            
end

[x,y] = meshgrid( x, y);
x = x(:); y = y(:);
        
            