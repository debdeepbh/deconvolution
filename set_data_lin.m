% generate theoretical data for multi-antenna model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Simulated data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_level = 5;

%% the impulse response to be used
%K = load('data/notches_260_0_0.txt');
%K = K(:,2)';

% get the theoretical data
testvec = testvec_gen(6);	% the spectrum of this one overlaps with that of K


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ANITA data                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a set of impulse responses
% borrow the impulse responses of the 5th signal
% drop ax, take aximp as the impulse responses
[ax, aximp] = prepsig(5);


%% Trim the observed signals to get 1024 length 


%%% %%%% taking all the impulse responses to be the same
%for i=1:15
%	aximp(i,:) = aximp(1,:);
%end


for i=1:15
	testconv(i,:) = conv(testvec, aximp(i,:))(1:1024);

	%%%% using the same noise level
	noiseax(i,:) = randn([1 length(testconv(i,:))])* noise_level;

	%%%%% noise of standard deviation 5i
	%noiseax(i,:) = randn([1 length(testconv(i,:))])*2*i;

	%if i<=3
	%	noiseax(i,:) = randn([1 length(testconv(i,:))])*3;
	%else
	%	noiseax(i,:) = randn([1 length(testconv(i,:))])*3*i;
	%end

	%%% Blurred noisy simulated data
	wax(i,:) =  testconv(i,:) + noiseax(i,:);
end


