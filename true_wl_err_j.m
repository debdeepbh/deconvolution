function [err_w_j, ratiounthres_j] = true_wl_err_j(j, scaling_j, faxobs, faximp, fori, w_ori, type, p, method, noiseax, sigma_leaked_j, rho) 
% given scaling_j, return the wavelet level error of the level j
% input sigmal(j) as sigmal_leaked_j

%% apply schiske
[fdec, mult] = fdecschiske(faxobs, faximp, fori, noiseax, scaling_j);


% get the original wavelet coeffs in the j level
w_ori_j = coeff(w_ori, p, j);

% compute the inverse wavelet transform of the schiske deconvolution
dec = real(ifft(fdec));
wdec = wtrans(dec, type, p);

% collect the coeffients of the j-th level
w_coeff = coeff(wdec, p, j);


% Applying a scaling on leaked noise variance
thrvec_j = sigma_leaked_j.*rho;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                    Level-dependent thresholding                     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch method
	case 'hard'
		mask = (abs(w_coeff) > thrvec_j);
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
	case 'oracle'
		w_ori_j  = coeff(w_ori, p, j);	% original signal
		mask = (abs(w_ori_j) > thrvec_j);
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
	case 'soft'
		mask = (abs(w_coeff) > thrvec_j);
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
		%% Soft thresholding
		coeff_thr = sign(coeff_thr) .* ( abs(coeff_thr) - thrvec_j);
	
	otherwise
		error('Wrong thresholding method');
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			    Level-dependent thresholding                     %
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	mask = (abs(w_coeff) > thrvec_j);
%	% Hard thresholding on the level
%	% If less than the noise level, drop it. Otherwise, keep it.
%	coeff_thr = w_coeff .* mask;
%
%	% Soft thresholding

	%% Number of un-thresholded coefficients
	ratiounthres_j = sum(mask)/length(mask); 

% compute the error of this wavelet level
err_w_j = norm(w_coeff - w_ori_j);
