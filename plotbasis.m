% plot the basis functions
function plotbasis(type, p, N)
B = getbasismat(type,p,N);
figure;
for j=1:p+1
	subplot(p+1,1,j);
	plot(0:N-1, shift(B(j,:),N/2));
	xlim([0 N-1]);
end

