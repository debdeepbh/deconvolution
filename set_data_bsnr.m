% generates simulated data for multi-antenna model
% fixed blurred-signal-to-noise ratio for each channel
function [testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_bsnr(testvec_num, bsnr) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 BSNR = 10 log_10 ( E|XH|^2/ E|N|^2)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bsnr/10 = log_10 ( E|XH|^2/ (N sigma^2))
% E|XH|^2/ (N sigma^2) = 10^(bsnr/10)
% N sigma^2 = E|XH|^2/(10^(bsnr/10))
% sigma^2 = E|XH|^2/ (N * 10^(bsnr/10))


% get the theoretical data
testvec = testvec_gen(testvec_num);	% the spectrum of this one overlaps with that of K
f_testvec = fft(testvec);

N = length(testvec);

for i=1:15
	aximp(i,:) = boxcar(i, N);

	testconv(i,:) = realconv(testvec, aximp(i,:));

	noise_level = sqrt( var( fft(testconv(i,:)) ) / (N * 10^(bsnr/10))); 

	%% Noise generation
	% random seed
	randn('seed', i);
	%%%% using the same noise level
	noiseax(i,:) = randn([1 length(testconv(i,:))])* noise_level;

	%%%%% noise of standard deviation 2i
	%noiseax(i,:) = randn([1 length(testconv(i,:))])*2*i;

	%%% Blurred noisy simulated data
	wax(i,:) =  testconv(i,:) + noiseax(i,:);

	%%%% fft of observed signal
	f_wax(i,:) = fft(wax(i,:));
	f_aximp(i,:) = fft(aximp(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           single channel                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick the best one, i.e. with highest signal-to-noise ratio
channel = 7;

testobs = wax(channel,:);
testimp = aximp(channel,:);
testnoise = noiseax(channel,:);

f_testobs = fft(testobs);
f_testimp = fft(testimp);

