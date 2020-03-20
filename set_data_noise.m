% generates simulated data for multi-antenna model
% fixed noise standard deviation for each channel
function [testvec, f_testvec, aximp, f_aximp, wax, f_wax, noiseax, testconv, testobs, testimp, testnoise, f_testobs, f_testimp] = set_data_noise(testvec_num, noise_level) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Simulated data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get the theoretical data
testvec = testvec_gen(testvec_num);	% the spectrum of this one overlaps with that of K
f_testvec = fft(testvec);

N = length(testvec);


for i=1:15
	aximp(i,:) = boxcar(i, N);
	%aximp(i,:) = boxcar(5, N);

	testconv(i,:) = realconv(testvec, aximp(i,:));


	%% Noise generation
	% random seed
	randn('seed', i);
	%%%% using the same noise level
	noiseax(i,:) = randn([1 N])* noise_level;
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
channel = 15; 

testobs = wax(channel,:);
testimp = aximp(channel,:);
testnoise = noiseax(channel,:);

f_testobs = fft(testobs);
f_testimp = fft(testimp);

