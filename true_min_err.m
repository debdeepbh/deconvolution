%function [w_thr, ratiounthres, thrvec]  = schiskeforwd_alt(faxobs, faximp, fori, B, type, p, noiseax, scaling, rho, method)
function [optsc_for_true_err, true_err] = true_min_err(faxobs, faximp, fori, B, type, p, noiseax, rho, method) 

% will store the argmin of true error
optsc_for_true_err = zeros(1, p+1);
true_err = zeros(1, p+1);

% prelim
ratiounthres = zeros(1, p+1); 
[M, N] = size(faxobs);

%% compute L and r
sqf_aximp = abs(faximp).^2;
L = sum(sqf_aximp, 1);
r = N*(noiseax^2)./(abs(fori).^2);	% using scalar noise sd to compute noise power

% wavelet transform of the original signal
ori = real(ifft(fori));
wori = wtrans(ori, type, p);

tic 

for j = 1:p+1

	% print j
	j

%%% find the scaling parameter that minimizes the true error
	sc_initial = 0.001;
	sc_step    = 0.01;
	sc_end     = 1;

	% brute force
	min_val = 10e4;
	opt_sc = 1;
	opt_err = 1;

	for sc_val=sc_initial:sc_step:sc_end

		%%% Get the leaked noise variance
		tempscaling = zeros(1,p+1);

		% put the current sc_val in the j-th place
		tempscaling(j) = sc_val;
		sigmal = getleakedsd(B, L, r, noiseax, tempscaling);
		% only take the j-th val
		sigma_leaked_j = sigmal(j);

	% compute the error
	y = true_wl_err_j(j, sc_val, faxobs, faximp, fori, wori, type, p, method, noiseax, sigma_leaked_j, rho);

		% compute the running minimum
		min_val = min(min_val, y);

		% if running minimum is the current val, store the scaling value
		if min_val == y
			opt_sc = sc_val; % update
			opt_err = y;
		end
	end

	%%% If the minimizer does not exist in [0,1]
	%%% we prefer to choose opt_sc =1 to be close to Schiske
	%if opt_sc == sc_step
	%	opt_sc = 1;
	%	opt_ratio = getratio_above_noise_schiske(1, B_j, L, r, noiseax, w_ori_j);

	%collect for each level j
	optsc_for_true_err(j) = opt_sc;
	true_err(j) = opt_err;




end

toc

