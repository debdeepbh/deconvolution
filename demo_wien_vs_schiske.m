%% Comparison between Wiener (for a single channel) and Schiske (for multi channel)
%% fixed noise sd used to create all channels
function [wien_err, schiske_err] = demo_wien_vs_schiske(testvec_num, noise_level, scaling) 

close all 

% set parameters
set_params

do_plot = 'no';

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

wien_est_all = zeros(size(wax));
wien_err = zeros(M,1);
for i=1:M
	%%% Single channel: run Wiener deconvolution
	f_testobs = f_wax(i,:);
	f_testimp = f_aximp(i,:);

	[fw, mult] = fdecwien(f_testobs, f_testimp, f_testvec, noise_level, scaling);
	wien_est = real(ifft(fw));

	wien_est_all(i,:) = wien_est;


	%wien_err(i) = norm(testvec - wien_est)/ norm(testvec);
	wien_err(i) = norm(testvec - wien_est);
end

%%% Multichanne: run Schiske deconvolution
% same noise level for each channel
[fschiske, mult_schiske] = fdecschiske(f_wax, f_aximp, f_testvec, noise_level, scaling);
schiske_est = real(ifft(fschiske));

%schiske_err = norm(testvec - schiske_est)/ norm(testvec);
schiske_err = norm(testvec - schiske_est);



switch do_plot
	case 'yes'
		figure
		plot(wax'); 
		xlim([1 N])

		figure
		plot(testvec);
		hold on
		plot(wien_est);
		%plot(schiske_est);
		xlim([1 N])
		hold off
		legend('Original', 'Estimate: Wiener', 'Estimate: Schiske')
	otherwise
		
end
