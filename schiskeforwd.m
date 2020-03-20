function [w, ratiounthres, thrvec]  = schiskeforwd(faxobs, faximp, fori, B, p, noiseax, scaling, rho, method)
% estimates the wavelet coefficients after doing deconvolution in Fourier then wavelet domain, based on supplied scaling values
% scaling should be of length p+1
% type is the wavelet filter type, p-th stage
% B is the wavelet basis constructed using getbasis()
% noiseax is scalar, constant noise sd across all channels

tic 

[M, N] = size(faxobs);

% if the supplied scaling is a scalar, vectorize is
if length(scaling) == 1
	scaling = zeros(1,p+1) + scaling;
end

%% compute L and r
sqf_aximp = abs(faximp).^2;
L = sum(sqf_aximp, 1);
r = N*(noiseax^2)./(abs(fori).^2);	% using scalar noise sd to compute noise power

% matlab allows empty matrix []
w = [];

for j=1:p+1
	Bsis = B(j,:);	% j-th row of basis matrix

	% do deconvolution using specific scaling for the j-th level
	%[fdec, mult] = fdecwien(fsig, fimp, fori, sigma, scaling(j));
	[fdec, mult] = fdecschiske(faxobs, faximp, fori, noiseax, scaling(j));


	%%% Modification of the index for the coarsest level
	%%% will use indexnow as a replacement of j
	indexnow = j;
	if j == p+1
		indexnow = p;
	end


	% prepare beta(l), the l-th wavelet coefficient of Fourier deconvolved signal
	beta = zeros(1,N/(2^indexnow));

	for l=0:(N/(2^indexnow) - 1)	% range of l, the localization index for current level
		% this is Psi_{j,l}, the (j,l)-th basis element
		Psi = shift(Bsis,(2^indexnow)*l);

		% dot product with the Fourier-based deconvolution, the wavelet coefficient
		% matlab dot product, non-commutative, u.v = sum{conj(u),v}
		% Dividing by N since dot(a, b) =  dot(fft(a), fft(b))/N 
		beta(l+1) = dot(fft(Psi), fdec)/N;
	end

	w = [w beta];
	%% Why is this complex number??
	w = real(w);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%                        Leaked noise variance                        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% computing the standard deviation of the leaked noise at the j-th (NOT indexnow) level
	% Dividing by N since dot(a,b) = dot(fft(a),fft(b))/N	% check on matlab
	%sigmal(j) = sqrt(sigma^2 * dot( (abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2)/N);

	if (length(noiseax) == 1)	% for scalar noise s.d.
		%sigmal(j) = sqrt(sigma^2 * dot( (abs(fft(Bsis))./abs(fimp)).^2, abs(mult).^2)/N);
		to_mult = L ./ ((L + r * scaling(j)).^2);
		%size(to_mult)	% check if a row vector
		dotprod_val = sum( to_mult .* (abs(Bsis).^2));	% sum over the freq
		%size(dotprod_val)	% check if a scalar
		sigmal(j) = sqrt( noiseax^2/N * dotprod_val);
	else
		%%% For vector valued sigma
		%sigmal(j) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Bsis))./abs(fimp)).^2, abs(mult).^2 )/(N.^2));
		error('Noise level is vector.')
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An alternate implementation will be to compute the Wiener deconvolution 
% for each j=1:p+1, compute wtrans and collect the j-th level and apply the
% threshold on it, but it involves computing irrelevant wavelet parts which 
% are thrown away in the final construction of w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% Plot %%%%%%%%%%%%%
% plot before thresholding
%plotcoeffs(w,p)
%print('before','-dpng')

% print before applying threshold
%figure;
%plot(iwtrans(w,type,p))
% plotthr(w,p,thrvec);
%title('before applying threshold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Level-dependent thresholding                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying a scaling on leaked noise variance
thrvec = sigmal.*rho;

% Applying the thresholding (in-place)
[w, ratiounthres, wnoise] = applythres(w, method, p, thrvec);
%w = keeplarge(w, 2);

toc

% plot after 
%plotcoeffs(w,p)

