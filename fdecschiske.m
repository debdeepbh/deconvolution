% Schiske deconvolution in the frequency domain
% Input: 	faxobs: Fourier transform of the noisy observed signal matrix
%		faximp: Fourier transform of the impulse response matrix
%		fori: Fourier transform of the guess for the original signal
% 			A possible guess will be fobs./sqrt(var(fimp))
%		noiseax: a row vector containing scalars for noise standard deviation, or the noise data matrix (time domain)
%		scaling: A regularized parameter to reduce the ringing, default =1
% Output:	fw:	Fourier transform of the estimate
%		mult:	Fourier shrinkage multiplier matrix

function [fw, mult] = fdecschiske(faxobs, faximp, fori, noiseax, scaling)


[M, N] = size(faxobs);

if ~exist('scaling')	% scaling is optional
	scaling = 1;
end

% assume that the fobs and fimp are of the same length

% right now, the noise value has to the same for each channel
[nsz_1, nsz_2] = size(noiseax);
if ( nsz_1 ~=1)
	error('The noise level has to be the same for all channels for now. Normalize the impulse responses to get a')
else
	if (nsz_2 ==1)	% scalar input
		sigmasq = N*noiseax^2;
	else			% it is the row vector 
		sigmasq = std(noiseax);
		%sigmasq = abs(fft(noise)).^2;
	end
end
	

%%% When we accommodate different noise levels for each channel
%if (nsz_2 == 1)	% either scalar or a column vector
	%if (nsz_1 == 1)	% noiseax is scalar
		%sigmasq = N*noiseax^2;
	%elseif		% noiseax is row vector, each row is the noise sd of a channel
		
	%end
%else	% it is the raw noise matrix
	%sigmasq = abs(fft(noise)).^2;
%end


% naive deconvolution in the fourier domain
% skip to avoid dividing by zero
%wfft = fobs./fimp;

sqf_aximp = abs(faximp).^2;

L = sum(sqf_aximp, 1);




r = sigmasq./(abs(fori).^2);

mult = sqf_aximp ./ (L + scaling * r);	% mult is a matrix, ith row is lambda_i

K = conj(faximp) ./ (L + scaling * r);	% K is a matrix

% estimate
fw = sum(faxobs .* K, 1);	% fw(f) = \sum_{i=1}^M Y_i(f) * K_i(f)
