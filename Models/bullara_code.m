
function L = bullara_code(h, lx, N, iter, Iband, flag)
	% Purpose:		simulates Bullara-DeDecker model
	% Notation:		S=precursor; X=xanthophore, M=melanophore; I=iridophore
	% Arguments:	h:    	width of M or X stripes
	% 				lx:   	long-range inhibition parameter
	%				[lx = 0.5 homogeneous; 1.5 black spots; 2.5 stripes; 3.5 yellow spots]
	%				N:    	lattice of size NxN
	%				iter: 	#iterations
	%				Iband:	width of centered horizontal iridophore band 
	%				flag: 	plots cells iff flag=1
	% Output:		NxN matrix L_ij in {0,1,2}={S,X,M} of final cell types on lattice
	% Example:		L = bullara_code(16, 2.5, 100, 1e7, 0, 1);

	%% Initialize random number generator for reproducibility
	%rng(1);

	%% Set transitions rates and interaction ranges
	par.bx = 1;			% birth rate of X
	par.bm = 0; 		% birth rate of M
	par.dx = 0; 		% death rate of X
	par.dm = 0; 		% death rate of M
	par.sx = 1; 		% inhibition rate of M by nearby X
	par.sm = 1; 		% inhibition rate of X by nearby M
	par.lx = lx; 		% rate  of long-range activation of M by X
	par.h  = h; 		% range of long-range activation of M by X

	%% Set lattice and simulation configuration
	par.N     = N;		% size of the grid
	par.iter  = iter;	% number of iterations for the simulation
	par.Iband = Iband;	% width of iridophore band

	%% Initialize lattice: S=0, X=1, M=2
	L = zeros(par.N,par.N);

	%% Initialize rates for S, X, and M transitions
	probs = [par.bx, par.bm, par.lx, 0, 0, 0, 0; ...	% transitions from S
			 0, 0, 0, par.dx, par.sm, 0, 0; ...			% transitions from X
			 0, 0, 0, 0, 0, par.dm, par.sx];			% transitions from M

	%% Initialize neighbors at distance 1
	nearest_neighbors = [1, 0; -1, 0; 0, 1; 0, -1];

	%% Random sample of cell locations, neighbors, and transitions
	cell = randi(par.N, par.iter, 2);
	neighbors = mod(cell + nearest_neighbors(randi(4,par.iter,1),:) - 1, par.N) + 1;
	trans = zeros(par.iter, 3);
	for j=1:3
		[~,trans(:,j)] = max(cumsum(probs(j,:))/sum(probs(j,:))>=rand(par.iter,1), [], 2);
	end

    %% Run simulation by iteratively updating transitions
	for t=1:par.iter
		switch trans(t, L(cell(t,1), cell(t,2)) + 1)
			case 1 	% S->X
				L(cell(t,1), cell(t,2)) = 1;
			case 2 	% S->M 			unless in iridophore band
				if (abs(cell(t,1)-par.N/2)>=par.Iband)
					L(cell(t,1), cell(t,2)) = 2;
				end
			case 3 	% S+Xh->M+Xh 	unless in iridophore band
				phi = 2*pi*rand(1,1);
				ind = round(cell(t,:) + par.h*[cos(phi),sin(phi)]);
				ind = mod(ind-1, par.N) + 1;
				if (L(ind(1), ind(2)) == 1 && abs(cell(t,1)-par.N/2)>=par.Iband)
					L(cell(t,1), cell(t,2)) = 2;
				end
			case 4 	% X->S
				L(cell(t,1), cell(t,2)) = 0;
			case 5 	% X+M1->S+M1 
				if (L(neighbors(t,1), neighbors(t,2)) == 2)
					L(cell(t,1), cell(t,2)) = 0;
				end
			case 6 	% M->S
				L(cell(t,1), cell(t,2)) = 0;
			case 7 	% M+X1->S+X1
				if (L(neighbors(t,1), neighbors(t,2)) == 1)
					L(cell(t,1), cell(t,2)) = 0;
				end
		end
	end

	%% Plot cells
	if (flag==1)	
		figure(1)
		imshow(L, [0,2], Colormap=[1 1 1; 1 1 0; 0 0 0], InitialMagnification=125e3/par.N)
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
