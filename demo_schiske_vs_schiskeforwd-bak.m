%% Comparison between Wiener (for a single channel) and Schiske (for multi channel)
%% Fixed bnsr is used to create the data
function [wien_err, schiske_err] = demo_wien_vs_schiske_bsnr(testvec_num, bsnr) 



% test original signal
%  1: two chirps
%  2: boxes
%  3: sawtooth/triangle
%  4: chirp
%  5: smooth impulsive
%  6: narrow impulsive

% default regularization parameter
scaling = 1;

% set parameters
set_params

%%% set the data

%[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_noise(testvec_num, noise_level);

[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_bsnr(testvec_num, bsnr);

% row vector with noise sd of channel
sd_noise = std(noiseax')';	


%%% Single channel: run Wiener deconvolution
[fw, mult] = fdecwien(f_testobs, f_testimp, f_testvec, testnoise, scaling);
wien_est = real(ifft(fw));

%%% Multichannel: run Schiske deconvolution
% different noise level for each channel
% normalize to minimize noise level to 1
f_wax_n = f_wax ./ sd_noise;
f_aximp_n = f_aximp ./ sd_noise;

[fschiske, mult_schiske] = fdecschiske(f_wax_n, f_aximp_n, f_testvec, 1, scaling);
schiske_est = real(ifft(fschiske));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                error                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wien_err = norm(testvec - wien_est)
schiske_err = norm(testvec - schiske_est)



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
plot(schiske_est);
xlim([1 1024])
hold off
legend('original', 'estimate: wiener', 'estimate: schiske')


