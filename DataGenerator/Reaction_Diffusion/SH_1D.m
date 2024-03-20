function u = SH_1D(mu, nu, viz)

% solves Swift-Hohenberg equation using ETD1 scheme
% computation is based on v = fft(u)
% u_t = -(Δ+1)^2 u - μu + νu^2 - u^3 on [-L,L] with periodic BCs

% nu = 1.6;
% mu = 0.2;

T = 50;

if nargin < 3
    viz = 0;
end

N  = 512;
L  = 32*pi;
dt = 0.1;
dx = 2*pi/N;
x  = dx*(1:N)';
x  = L*(x-pi)/pi;

%% initial condition
L1 = 4*pi;
u = (-tanh(x-L1) + tanh(x+L1)).*cos(x) + 0.0*rand(1,N)';

%% compute linear operator
k   = [0:N/2 -N/2+1:-1]*(pi/L);           % wave numbers
k = k';
LU  = -(-k.^2+1).^2 - mu;  				% linear operator
EXP = exp(LU*dt);                       % Exact linear bit
ETD = (exp(LU*dt)-1)./LU;               % ETD1 coeffs of linear operator

uT = fft(u);                           % Time step done in Fourier space


for t=0:dt:T                           % evolve via ETD  method
	f = nu*u.^2 - u.^3;                 % nonlinear RHS
    
	% exponential timestep in fourier space
	fT = fft(f);
	uT = EXP.*uT + ETD.*fT;
	u  = real(ifft(uT));
    
	if max(max(abs(u)))>100; return; end;
    if viz
        figure(1)
        plot(x,u); ylabel('u'); xlabel('x'); drawnow;
    end
    
    
end















