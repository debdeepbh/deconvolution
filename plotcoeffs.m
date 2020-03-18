% plot the coefficients of wavelet transform w, the p-th stage wavelet transform
function plotcoeffs(w,p)
N = length(w);
type = '.';

figure;

[ha, pos] = tight_subplot(p+1,1,[.01 .03],[0.05 .01],[.1 .05]);

for i=1:p
	%subplot(p+1,1,i);
	axes(ha(i));
	gap = 2^i;
	vals = coeff(w,p,i);
	stem(1:gap:N, vals,type)
	xlim([1 N]);
	ylabel(strcat('Level: ',num2str(i)));

	%axis([1 N -max(abs(vals)) max(abs(vals))]);
	% to remove the x-axis labels, markers
end

%plot the last coarsest part, with the same gap as last
axes(ha(p+1));
%subplot(p+1,1,p+1);
vals = coeff(w,p,p+1);
stem(1:gap:N, vals, type)
xlim([1 N]);
ylabel(strcat('Coarsest level:',num2str(p)));
%axis([1 N -max(abs(vals)) max(abs(vals))])

% remove x-label except for the last one
set(ha(1:p),'XTickLabel',''); 
%set(ha,'YTickLabel','');
