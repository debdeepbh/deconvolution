% pad frequencies to upsample to get a vector of length N
% can also be used to window without padding, when N = length(w)
% window then pad to avoid ringing,
% eps is the proportion of *actual* frequency of the original signal to get modified
% eps=1 implies window all, ie planktaper eps = 0.5
function z = padfreq(w, N, eps)
% caution: length of w must be even
n = length(w);
if (n >N)
	disp('nothing to pad here');
	z = w;
else
	fw = fft(w);
	extral = N-n;

	% a plank-taper window
	ptw = planktaper(n, eps/2);
	cutoff = [ptw(floor(n/2)+1:n) ptw(1:floor(n/2))]; % cut in half and flip
	
	% window original fft
	newf = fw.*cutoff;

	% zero-padding
	newf = [newf(1:floor(n/2)) zeros(1,extral) newf(floor(n/2)+1:n)];

	z =real(ifft(newf));
end



