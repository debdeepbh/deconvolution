function [true_scaling_min, min_err_all] = true_min_err_alt(faxobs, faximp, fori, B, type, p, noiseax, rho, method) 


% wavelet transform of the original signal
ori = real(ifft(fori));
wori = wtrans(ori, type, p);


tic


%% define error function, relative error
err_fun = @(scaling) norm(iwtrans(schiskeforwd_alt(faxobs, faximp, fori, B, type, p, noiseax, scaling, rho, method), type, p) - ori)/norm(ori);


%%%% fminsearch

 scaling_guess = zeros(1, p+1) + 0.5;
 %scaling_guess = zeros(1, p+1) + 0.01;

%scaling_guess = [0.024391   0.024391   0.149141   0.663734];

tic

[true_scaling_min, min_err_all]  = fminsearch(err_fun, scaling_guess);

toc
