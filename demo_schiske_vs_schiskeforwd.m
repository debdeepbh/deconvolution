%% Comparison between Schiske and schiskeforwar (for multi channel)
%% fixed noise sd used to create all channels
function [schiske_err, schiskeforwd_err] = demo_schiske_vs_schiskeforwd(testvec_num, noise_level, scaling) 

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

%noisetype = 'constant';
noisetype = 'bsnr';

% rho
rho = 1;

% Thresholding method
method = 'hard';

% set parameters
set_params

%%% get the data

switch noisetype
	case 'constant'
		noisetype

		[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_noise(testvec_num, noise_level);
	case 'bsnr'
		noisetype

		bsnr = noise_level;
		[testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_bsnr(testvec_num, bsnr);

	% row vector with noise sd of channel
	sd_noise = std(noiseax')';	

		% different noise level for each channel
		% normalize to minimize noise level to 1

		wax     = wax    ./ sd_noise;
		aximp   = aximp  ./ sd_noise;
		f_wax   = f_wax  ./ sd_noise;
		f_aximp = f_aximp./ sd_noise;

		% replace the noise_level in the end
		noise_level = 1;
	otherwise
		error('Incorrect noise type');
		
end



% Single channel: run Wiener deconvolution
%[fw, mult] = fdecwien(f_testobs, f_testimp, f_testvec, noise_level, scaling);
%wien_est = real(ifft(fw));

%%% Multichannel: run Schiske deconvolution
% same noise level for each channel
% Take scaling = 1 here to do a true schiske
[fschiske, mult_schiske] = fdecschiske(f_wax, f_aximp, f_testvec, noise_level, 1);
schiske_est = real(ifft(fschiske));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% True error producing scaling parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[optsc_for_true_err, true_err] = true_min_err(f_wax, f_aximp, f_testvec, B, type, p, noise_level, rho, method);
%
%optsc_for_true_err
%true_err
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Multichanel: run schiskeforwd deconvolution
% old implementation
%[w, ratiounthres, thrvec]  = schiskeforwd(f_wax, f_aximp, f_testvec, B, p, noise_level, scaling, rho, method);

% alternate implementation
[w, ratiounthres, thrvec]  = schiskeforwd_alt(f_wax, f_aximp, f_testvec, B, type, p, noise_level, scaling, rho, method);
schiskeforwd_est = iwtrans(w, type, p);

ratiounthres
thrvec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                error                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
schiske_err = norm(testvec - schiske_est)
schiskeforwd_err = norm(testvec - schiskeforwd_est)

perc_improvement = (schiske_err - schiskeforwd_err)/ schiske_err * 100




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Plot                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot all noisy blurred signals
%figure
%plot(wax'); 
%xlim([1 N])

figure
plot(testvec);
hold on
plot(schiske_est);
plot(schiskeforwd_est);
xlim([1 N])
hold off
legend('original', 'estimate: schiske', 'estimate: schiskeforwd')


