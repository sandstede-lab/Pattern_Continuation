function [u,x,y] = SH_2D(mu, nu, Uinit, viz, N, T)
%GENERATE_SH_TIMELAPSE Summary of this function goes here
%   Displays and saves timelapses of the Swift-Hohenberg model at specified timepoints
% solves Swift-Hohenberg equation using ETD1 scheme
% computation is based on v = fft(u)
% u_t = -(Δ+1)^2 u + μu + νu^2 - u^3 on [-L,L] with periodic BCs


% hexagons are stable for -nu^2/15 < mu < 16*nu^2/3
% rolls are stable for 4*nu^2/3 < mu

%nu = 0.5;			% equation parameter: beta    in my notes
%mu = 2*nu^2/3;		% equation parameter: epsilon in my notes
%mu = nu^2/3;      % bifurc of mixed modes/rolls - appears like hexagons
%                       this one is actually ok around small mu, but the
%                       smaller mu gets, the longer it takes to integrate
%mu = 4 * nu^2/3;  % bifurc of mixed modes/false hexagons - appears labyrinthine
%mu = - nu^2/60;   % bifurcation of hexagons - appears like hexagons stabilize slowly from one end to another
%nu = 0.5;
%mu = 4 * nu^2/3;
% (0.01, 0.01), T = 3000 -> stripes
% reference_pt = [0.1230, 0.3230];
% unit_delta = [0.2144, -0.9767];
% pt = reference_pt + scale * unit_delta;
% mu = pt(1);
% nu = pt(2);
% (0.01, sqrt(3 * mu)), T = 3000 -> hexagons
% (0.01, sqrt(3 * mu)) - sqrt(mu)), T = 3000 -> stripes sometimes
% (0.008, sqrt(3 * mu)) - sqrt(mu)), T = 3000 -> stripes sometimes
%nu = sqrt(3 * mu / 4);

%P_h0 = -8*s^2/135;               % limit eps between hexagons and 0 stable solutions
%P_hr = (5*s^2/15)*(7+3*sqrt(6)); % limit eps between hexagons and rolls stable solutions



if nargin < 4
    viz = 0;
end
if nargin < 6
    
    T = 1000;%/2*1.5;
end
% T = 600;
if nargin < 5
    N  = 128;
end
L  = 16*pi;
dt = 0.1;
dx = 2*pi/N;
x  = dx*(1:N)';
x  = L*(x-pi)/pi;
y  = x';
[xx,yy] = meshgrid(x,y);

%% initial condition
if nargin < 3
    u = 1.5*rand(size(xx));
else
    switch Uinit
        case 'spots'
            u = 0.5/3*(cos(xx) + cos((xx+sqrt(3)*yy)/2) + cos((xx-sqrt(3)*yy)/2));
        case 'stripes'
            u = 0.2*cos(xx);
        case 'random'
            u = 1.5*rand(size(xx));
    end
end

% if init == 1
% %% hexagons
cc  = -20;
% u = 2/3*(cos(xx-cc) + cos(((xx-cc)+sqrt(3)*(yy-0))/2) + cos(((xx-cc)-sqrt(3)*(yy-0))/2));
% u = 0.5/3*(cos(xx) + cos((xx+sqrt(3)*yy)/2) + cos((xx-sqrt(3)*yy)/2));
% u = u + 0.1*rand(size(xx));
% 
% %% rolls
% u(65:128, :) = 0.2*cos(xx(65:128, :));
% end
% 
% if init == 0
% %% random
% u = 1.5*rand(size(xx));
% end

% compute linear operator
k = [0:N/2 -N/2+1:-1]*(pi/L);           % wave numbers
[kkx,kky]=meshgrid(k,k);
LU = - (-(kkx.^2 + kky.^2)+1).^2 + mu;  % linear operator
EXP = exp(LU*dt);                       % Exact linear bit
ETD = (exp(LU*dt)-1)./LU;               % ETD1 coeffs of linear operator

uT = fft2(u);                           % Time step done in Fourier space

j=0;
t_show = round(T*linspace(0,1,9) );

for t=0:dt:T    
%      t% evolve via ETD  method
     f = nu*u.^2 - u.^3;                 % nonlinear RHS

     % exponential timestep in fourier space
     fT  = fft2(f); uT = uT.*EXP + fT.*ETD;      
     u   = real(ifft2(uT));

	 if max(max(abs(u)))>100; return; end
     if viz
         if t == t_show(1) || t == t_show(2) || t == t_show(3) || t == t_show(4) || t == t_show(5) || t == t_show(6) || t == t_show(7) || t == t_show(8) || t == t_show(9)
            j = j + 1;
            surf(xx,yy,u); ylabel('y'); xlabel('x'); view([0 90]); title(['t=' num2str(t_show(j)) ' max u=' num2str(max(max(u)))]);
            shading interp; colorbar;
            pause(0.1);
    %         saveas(gcf, [save_location '/scale=' num2str(abs(scale)) '/m=' num2str(mu) '_n=' num2str(nu) '_t=' num2str(t) '.png']);
         end
     end

%      if mod(t,disp_interval) == 0                 % display at intervals
%         title(['t=' num2str(t) ' max u=' num2str(max(max(u)))]);
% 		subplot(1,2,1);
%         surf(xx,yy,u); ylabel('y'); xlabel('x'); view([0 90]);
%         shading interp; colorbar; drawnow;
% 		subplot(1,2,2);
%         surf(xx,yy,v);view([-90 90]);
%         shading interp; axis off; colorbar; drawnow;
%      end
     
     % sample from u surface for taking
     % see mu substitutions above - use mu and nu as variable parameters
     % when doing bifurcation tracing
end


x = xx(:);
y = yy(:);

end