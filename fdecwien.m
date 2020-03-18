% Wiener deconvolution in the frequency domain
% Input: 	fobs: Fourier transform of the noisy observed signal
%		fimp: Fourier transform of the impulse response	
%		fori: Fourier transform of the guess for the original signal
% 			A possible guess will be fobs./sqrt(var(fimp))
%		noise: a scalar for noise standard deviation, or the noise data (time domain)
%		scaling: A regularized parameter to reduce the ringing, default =1
% Output:	fw:	Fourier transform of the estimate
%		mult:	Fourier shrinkage multiplier

function [fw, mult] = fdecwien(fobs, fimp, fori, noise, scaling)


N = length(fobs);
if ~exist('scaling')	% scaling is optional
	scaling = 1;
end

% assume that the fobs and fimp are of the same length

% naive deconvolution in the fourier domain
% skip to avoid dividing by zero
%wfft = fobs./fimp;

	% definition
	hsq = abs(fimp).^2;

	% construct the noise power
	if (length(noise)==1)	% noise is actually the variance
		sigmasq = N*noise^2;
	else			% it is the raw noise
		sigmasq = abs(fft(noise)).^2;
	end

	%mult = hsq ./( hsq + scaling*N*(sigma.^2)./(abs(fori).^2));
	%mult = hsq ./( hsq + scaling*(sigmasq)./(abs(fori).^2));
	%mult = hsq ./( hsq + scaling*(sigmasq)./((abs(fori).^2) .*(hsq) ));

	fw = fobs .* conj(fimp) ./( hsq + scaling*(sigmasq)./(abs(fori).^2));

	% to correct the error of considering y as x, we multiply by the norm of h  so that ||xtil|| = ||x|| 
	% where xtilhat = yhat/hhat*||hhat||, hat is Fourier transform
	% This gives a much better result while doing wiener
	%fw = fobs .* conj(fimp) ./( hsq + scaling*(sigmasq)./((abs(fori).^2)./var(fimp)) );

	% we need to output multi as well, for wienforwd.m
	mult = hsq ./( hsq + scaling*(sigmasq)./(abs(fori).^2) );
