% given the p-th stage wavelet transform, get the qth level coeff
% q can be 1:p+1. For q=p+1, output the coarsest part
% Slow implementation: generate dyadic indices and call directly instead
function values = coeff(w, p, q)

maxet = q;
%if q=p+1, print the coarsest part
if q == p+1
	% so iterate p times only
	maxet = p;
end

% a recursive vector
recv = w;
for i=1:maxet
	N = length(recv);
	% first half of the vector
	first = recv(1:N/2);
	% second half, discard it
	sec = recv(N/2+1:N);
	% set the second part as the recursive vector
	recv = sec;
end

% it is the first part that we want
values = first;

%if q=p+1, print the coarsest part
if q == p+1
	values = sec;
end



