% generate theoretical data for multi-antenna model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Simulated data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_level = 2;

% get the theoretical data
testvec = testvec_gen(2);	% the spectrum of this one overlaps with that of K
f_testvec = fft(testvec);

aximp = zeros(15, 1024);

[ax, aximp_anita] = prepsig(5);


for i=1:15
	aximp(i,:) = boxcar(i, 1024);

	testconv(i,:) = realconv(testvec, aximp(i,:));

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
channel = 15;
testobs = wax(channel,:);
testimp = aximp(channel,:);
testnoise = noiseax(channel,:);

f_testobs = fft(testobs);
f_testimp = fft(testimp);

