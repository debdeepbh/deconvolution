function [w_decomp] = decompose_wl(w, B) 

[p_p_1, N] = size(B);
p  = p_p_1 -1;

w_decomp = zeros(size(B));

for j = 1:p+1

	Bsis = B(j,:);

	% wavelet coefficients of current level
	w_coeff = coeff(w, p, j);

	if j==(p+1)
		indexnow = p;
	else
		indexnow = j;
	end
	
	decomp_j = zeros(1,N);

	for l=0:(N/(2^indexnow) - 1)	% range of l, the localization index for current level
		% this is Psi_{j,l}, the (j,l)-th basis element
		Psi = shift(Bsis,(2^indexnow)*l);

		decomp_j = decomp_j + Psi * w_coeff(l+1);

	end

	w_decomp(j,:) = decomp_j;
end


%% plot
[ha, pos] = tight_subplot(p+1,1, [.01 .03], [0.05 .01], [.1 .05]);
pos
ha
for j=1:p+1
	axes(ha(j));	% change plot number
	%tight_subplot(p+1, 1, j);
	plot(w_decomp(j,:));
	%axis off;
        %set(gca,'xtick',[])
    %set(gca,'xticklabel',[])
	xlim([1 N]);
	max_abs = max(abs(w_decomp(j,:)));
	ylim([-max_abs, max_abs]);
end

%% Check if adds up
%figure
%plot(sum(w_decomp))
