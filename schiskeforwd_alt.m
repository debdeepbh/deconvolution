function [w_thr, ratiounthres, thrvec]  = schiskeforwd_alt(faxobs, faximp, fori, B, type, p, noiseax, scaling, rho, method)
% alternative implementation of schiskeforwd.m
% Note: type is required here
% estimates the wavelet coefficients after doing deconvolution in Fourier then wavelet domain, based on supplied scaling values
% scaling should be of length p+1. If scaling==(-1) compute the optimum scaling parameter.
% type is the wavelet filter type, p-th stage
% B is the wavelet basis constructed using getbasis()
% noiseax is scalar, constant noise sd across all channels

do_plot = 'no';
%do_plot = 'yes';

ratiounthres = zeros(1, p+1); 

[M, N] = size(faxobs);

%% compute L and r
sqf_aximp = abs(faximp).^2;
L = sum(sqf_aximp, 1);
r = N*(noiseax^2)./(abs(fori).^2);	% using scalar noise sd to compute noise power


% wavelet transform of the original signal
ori = real(ifft(fori));
w_ori = wtrans(ori, type, p);

% if the supplied scaling is a scalar, vectorize is
if length(scaling) == 1
	if (scaling >= 0)
		scaling = zeros(1,p+1) + scaling;
	elseif (scaling == -1)	% then compute the optimum scaling
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%                  Get the optimum scaling a priori                   %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% rewriting the input scaling
		ratio_above = zeros(size(scaling));
		for j = 1:p+1
			B_j = B(j,:);	% basis element
			w_ori_j = coeff(w_ori, p, j);	% wavelet coeffs

			% optimal scaling
			[scaling(j), ratio_above(j)] = getopt_j(B_j, L, r, noiseax, w_ori_j, 'bisec');
		end
		printf('\n')
		scaling
		ratio_above
	elseif (scaling == -2)	% compute the minimizer of true MSE
		[true_scaling_min, min_err_true_mse] = true_min_err_alt(faxobs, faximp, fori, B, type, p, noiseax, rho, method)
		
		% stop execution
		return
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Get the leaked noise variance
sigmal = getleakedsd(B, L, r, noiseax, scaling);

% Applying a scaling on leaked noise variance
thrvec = sigmal.*rho;

% empty wavelet coeffs
w = [];
w_thr = [];

%%% Perform schiske deconv using scaling values
for j=1:p+1
	% do deconvolution using specific scaling for the j-th level
	[fdec, mult] = fdecschiske(faxobs, faximp, fori, noiseax, scaling(j));

	% compute the wavelet transform of the schiske deconvolution
	dec = real(ifft(fdec));
	wdec = wtrans(dec, type, p);

	switch do_plot
		case 'yes'
			plot(dec);
			hold on
		otherwise
			
	end

	% collect the coeffients of the j-th level
	w_coeff = coeff(wdec, p, j);
	% put it in the j-th level
	w = [w w_coeff];
	%printf('Schiske error in level %d with %f: %f\n', j, scaling(j), norm(dec - ori));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                    Level-dependent thresholding                     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch method
	case 'hard'
		mask = (abs(w_coeff) > thrvec(j));
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
	case 'oracle'
		w_ori_j  = coeff(w_ori, p, j);	% original signal
		mask = (abs(w_ori_j) > thrvec(j));
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
	case 'soft'
		mask = (abs(w_coeff) > thrvec(j));
		%% Hard thresholding on the level
		%% If less than the noise level, drop it. Otherwise, keep it.
		coeff_thr = w_coeff .* mask;
		%% Soft thresholding
		coeff_thr = sign(coeff_thr) .* ( abs(coeff_thr) - thrvec(j));
	
	otherwise
		error('Wrong thresholding method');
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% place it
	w_thr = [w_thr coeff_thr];

	%% Number of un-thresholded coefficients
	ratiounthres(j) = sum(mask)/length(mask); 

end
	switch do_plot
		case 'yes'
			hold off
			legend('1', '2', '3', '3+1');
		otherwise
			
	end

switch do_plot
	case 'yes'
		
		%% print before denoising
		%printf('Error before denoising: %f \n', norm( iwtrans(w, type, p) - ori));
		%

		plotthr(w_ori, p, thrvec)
		plotthr(w_thr, p, thrvec)
		%plotthr(w_thr, p, thrvec)

	otherwise
		
end


% Applying the thresholding (in-place)
%[w, ratiounthres, wnoise] = applythres(w, method, p, thrvec);
%[w, ratiounthres, wnoise] = threshold(w, method, p, thrvec);
%w = keeplarge(w, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Print                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigmal
%scaling


