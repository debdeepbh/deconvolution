% plot all the bases of p-stage wavelet transform
% outputs an p+1 by N matrix B
% columns of B are f_1, ..., f_p, g_p
% translates of which gives you the basis elements
% Psi_{-j,k} = R_{2^j * k} f_j, for j = 1,...,p 
% Phi_{-j,k} = R_{2^j * k} g_j, for j = p
% k = 0, ... , (N/(2^j)-1)
function B = getbasismat(type, p, N)
[u, v] = filt(type, N);

% empty vectors to store the intermediate values of u and v
U = V = zeros(p,N);

U(1,:) = u;
V(1,:) = v;
for k=2:p
	% fold it k times
	u_folded = u;
	v_folded = v;
	for j=1:k-1
		u_folded = fold(u_folded);
		v_folded = fold(v_folded);
	end
	% upsample it k times
	u_upped = u_folded;
	v_upped = v_folded;
	for j=1:k-1
		u_upped = up(u_upped);
		v_upped = up(v_upped);
	end
	% store the vectors
	U(k,:) = u_upped;
	V(k,:) = v_upped;
end

% empty vectors to store shifted basis elements
Psi = zeros(p, N);
Phi = zeros(p);


f_old = v;
g_old = u;

Psi(1,:) = f_old;

for j=2:p
	% f and g for this stage
	f_new = realconv(g_old, V(j,:));
	g_new = realconv(g_old, U(j,:));

	% psi_{-j,k} = R_{2^j.k} f_j
	% phi_{-j,k} = R_{2^j.k} g_j
	% shift the basis to make it centered at N//2
	% so that it looks nice while plotting
	% note that this is the particular basis element 
	% where 2^j . k = N/2, solve for k, j given
	Psi(j,:) = f_new;

	% loop condition
	f_old = f_new;
	g_old = g_new;
end

% the last Phi
Phi = g_old;

%%plotting the bases
%figure;
%for j=1:p
%	subplot(p+1,1,j); 
%	plot(0:N-1, shift(Psi(j,:),N/2));
%	axis([0 N-1]);
%end
%
%Phi = shift(g_old, N/2);
%subplot(p+1,1,p+1);
%plot(0:N-1, Phi);
%axis([0 N-1]);

B = [Psi; Phi];
