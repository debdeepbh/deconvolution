function [w, ratiounthr, wnoise] = applythres(wo, rule, p, thr)
% apply threshold on the p-th resolution wavelet transform
% thr should be a vector of length p+1
% wnoise is the dropped part of the original signal

% stores the proportion of wavelet coeffs that have been thresholded
ratiounthr = zeros(1,p+1);

% don't want to replace the original
w = wo;
N = length(w);

% check
if length(thr) ~= p+1
	disp('threshold vector length incorrect');
end


% staring index
starting = 1;
len = N/2;
for k=1:p
	% reset the counter of unthresholded coefficients
	thresholded = 0;

	ending = starting+len-1;

	for l=starting:ending

		switch rule
		case 'hard'
			if( abs(w(l)) < thr(k))
				w(l) = 0;
				% update ratio counter
				thresholded = thresholded + 1;
			end
		case 'soft'
			if( abs(w(l)) >= thr(k))
				w(l) = sign(w(l)) * (abs(w(l)) - thr(k));
			else
				w(l) = 0;
				% update ratio counter
				thresholded = thresholded + 1;
			end
		case 'wien'	% wavelet domain wiener shrinkage 
			% applies to all the wavelet coefficients
			% irrespective of the threshold
			% Here, the threshold should be the 
			% leaked noise variance at that level
			w(l) = (w(l)^2)/(w(l)^2 + thr(k)^2);
		otherwise
			disp('Wrong thresholding method.');
		end

	end

	% get the ratio of coefficients that were bigger than the 
	% noise variance
	Nj = length(starting:ending);
	ratiounthr(k) = (Nj - thresholded)/Nj;


	% for the next loop
	len = len/2;
	starting = ending + 1;
end

% for the coarsest level Phi_{-p,k}
% reset the counter of unthresholded coefficients
thresholded = 0;
for l=starting:N
	%%if( w(l) < thr(p+1))
	%%	w(l) = 0;
	%%	% update ratio counter
	%%	thresholded = thresholded + 1;
	%%end
	switch rule
	case 'hard'
		if( abs(w(l)) < thr(p+1))
			w(l) = 0;
			% update ratio counter
			thresholded = thresholded + 1;
		end
	case 'soft'
		if( abs(w(l)) >= thr(p+1))
			w(l) = sign(w(l)) * (abs(w(l)) - thr(p+1));
		else
			w(l) = 0;
			% update ratio counter
			thresholded = thresholded + 1;
		end
	case 'wien'	% wavelet domain wiener shrinkage 
		% applies to all the wavelet coefficients
		% irrespective of the threshold
		% Here, the threshold should be the 
		% leaked noise variance at that level
		w(l) = (w(l)^2)/(w(l)^2 + thr(p+1)^2);
	otherwise
		disp('Wrong thresholding method.');
	end
end

% get the ratio
Nj = length(starting:N);
ratiounthr(p+1) = (Nj - thresholded)/Nj;

% print the threshold values used
wnoise = wo - w;
