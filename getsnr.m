% get the snr of a signal based on
% peak-to-peak/(2*rms)
function [snrval, M, m, a, b] = getsnr(z)
%snrval = (max(z) - min(z))/(2* sum(abs(z).^2)/sqrt(length(z)));

% need to exclude peak-to-peak region
n = 10;	% how much around the peaks to exclude
M = max(z);
m = min(z);
maxind = find(z == M)(1);
minind = find(z == m)(1);
first = min(minind,maxind);
last = max(minind,maxind);
peakarea = last - first;
a = first - n*peakarea;
b = last + n*peakarea;
if a<1
	a =1;
end
if b>length(z)
	b = length(z);
end
restvec = [z(1:a) z(b:length(z))];

snrval = (M - m)/(2*sqrt(sum(abs(restvec).^2)/length(restvec)));
