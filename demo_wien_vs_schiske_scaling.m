%% Comparison between Wiener (for a single channel) and Schiske (for multi channel)
%% fixed noise sd used to create all channels
function [wien_err, schiske_err] = demo_wien_vs_schiske_scaling(testvec_num, noise_level, scaling) 

close all

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
[fw, mult] = fdecwien(f_testobs, f_testimp, f_testvec, noise_level, scaling);
wien_est = real(ifft(fw));

%%% Multichanne: run Schiske deconvolution
% same noise level for each channel
[fschiske, mult_schiske] = fdecschiske(f_wax, f_aximp, f_testvec, noise_level, scaling);
schiske_est = real(ifft(fschiske));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                error                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wien_err = norm(testvec - wien_est);
schiske_err = norm(testvec - schiske_est);



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

%% Noise seed consistency test
%disp('sum')
%sum(sum(noiseax))
