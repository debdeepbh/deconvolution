%% Comparison between Wiener (for a single channel) and Schiske (for multi channel)
%% fixed noise sd used to create all channels
function [wien_err, wienforwd_err] = demo_wien_vs_wienforwd_noise(testvec_num, noise_level, scaling) 





% test original signal
%  1: two chirps
%  2: boxes
%  3: sawtooth/triangle
%  4: chirp
%  5: smooth impulsive
%  6: narrow impulsive

% default regularization parameter
%scaling = 1;

% set parameters
set_params

%%% get the data

[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_noise(testvec_num, noise_level);

%[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_bsnr(testvec_num, bsnr);


%%% Single channel: run Wiener deconvolution
% Take scaling =1 here to do a true Wiener
[fw, mult] = fdecwien(f_testobs, f_testimp, f_testvec, noise_level, 1);
wien_est = real(ifft(fw));

%%% Single channel: run wienforwd deconvolution
[w, ratiounthres, thrvec]  = wienforwd(f_testobs, f_testimp, f_testvec, B, p, noise_level, scaling, rho, method);
wienforwd_est = iwtrans(w, type, p);

% print the unthrsholded ratio
ratiounthres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                error                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wien_err = norm(testvec - wien_est)
wienforwd_err = norm(testvec - wienforwd_est)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Plot                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(wax'); 
xlim([1 1024])

figure
plot(testvec);
hold on
plot(wien_est);
plot(wienforwd_est);
xlim([1 1024])
hold off
legend('original', 'estimate: wiener', 'estimate: wienforwd')


