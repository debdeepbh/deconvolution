function ratio_j = getratio_above_noise_schiske(scaling_j, B_j, L, r, noiseax, w_ori_j) 
% Given one scalar component of the scaling vector (and other info), compute the ratio of  coefficients above the leaked noise level
% w_ori_j is the j-th level wavelet coefficients of the original signal
% B_j is the j-th row of the basis matrix

N = length(B_j);

% get the leaked noise variances, given the scaling parameters
to_mult = L ./ ((L + r * scaling_j).^2);
dotprod_val = sum( to_mult .* (abs(fft(B_j)).^2) );	% sum over the freq
leaked_noise_level_j = sqrt( noiseax^2/N * dotprod_val);

% for a certain level, collect the ratio of wx coeffs that are bigger than the leaked noise
num_coef_above = sum (abs(w_ori_j) > leaked_noise_level_j);
ratio_j = num_coef_above/(length(w_ori_j));


