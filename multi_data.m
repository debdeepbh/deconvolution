% generate theoretical data for multi-antenna model

% the impulse response to be used
K = load('data/notches_260_0_0.txt');
K = K(:,2)';

% get the theoretical data
testvec_6;	% the spectrum of this one overlaps with that of K
testyori = [testvec zeros(1,512)];

% get a set of impulse responses
% borrow the impulse responses of the 5th signal
% drop ax, take aximp as the impulse responses
[ax, aximp] = prepsig(5);



% pad by zero to make longer
testyori = [testvec zeros(1,512)];

%%% %%%% taking all the impulse responses to be the same
%for i=1:15
%	aximp(i,:) = aximp(1,:);
%end


for i=1:15
	testconv(i,:) = conv(testvec, aximp(i,:))(1:1024);
	% noise of standard deviation 5i
	noiseax(i,:) = randn([1 length(testconv(i,:))])*2*i;
%% %%using same noise level
	%noiseax(i,:) = randn([1 length(testconv(i,:))])*11;
	%if i<=3
	%	noiseax(i,:) = randn([1 length(testconv(i,:))])*3;
	%else
	%	noiseax(i,:) = randn([1 length(testconv(i,:))])*3*i;
	%end
	wax(i,:) =  testconv(i,:) + noiseax(i,:);
end


% finally, we get 
% the original data testyori
%  wax, each row is an observed signal
% aximp, each row is an impulse response
%  noiseax, each row is the noise data
% 

