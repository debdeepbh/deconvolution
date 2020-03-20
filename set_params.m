%% set up the parameters to be used 
	 N = 512;
	 M = 15;
	 type='d10';
	 p=3;
	 method='hard';
	 B = getbasismat(type, p, N);
	 fB = fft(B')';

% running the same script again with one of the global variables changed does not change the value of the variable (if called without global)
	%global M = 15;
	%global type='d10';
	%%global type='meyer';
	%global p=3;
	%global rho=1;
	%global method='hard';
	%global B = getbasismat(type, p, N);
	%global fB = fft(B')';

	


