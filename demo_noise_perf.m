function demo_noise_perf(testvec_num) 
% noise performance

scaling = 1;

M = 15;

compare = 'schiske_vs_schiskeforwd';
%compare = 'wien_vs_schiske';

levels = 1:30;

s_e = zeros(size(levels));
sf_e = zeros(size(levels));

wien_err_all = zeros(M, length(levels));
schiske_err_all = zeros(size(levels));	% 1x lenght(levels)

ind=1;
for level=levels
	switch compare
		case 'schiske_vs_schiskeforwd'
			[s_e(ind), sf_e(ind)] = demo_schiske_vs_schiskeforwd(testvec_num, level, scaling);
		case 'wien_vs_schiske'


			[wien_err, schiske_err] = demo_wien_vs_schiske(testvec_num, level, scaling) 
			wien_err_all(:,ind) = wien_err;
			schiske_err_all(ind) = schiske_err;
		otherwise
			error('Wrong keyword for comparison')
			
	end
	ind=ind+1;
end

switch compare
	case 'schiske_vs_schiskeforwd'
		plot(levels, s_e, levels, sf_e)
		legend('Schiske', 'Multichannel ForWaRD')
		xlabel('BSNR (dB)')
		ylabel('Relative l^2 error')
	case 'wien_vs_schiske'
		for i=1:M
			plot(levels, wien_err_all(i,:))
			hold on

			% string for ith legend
			legendStr{i} = strcat('channel: ',num2str(i));

		end
		plot(levels, schiske_err_all, 'o-'); 
			%legend('Schiske', 'Multichannel ForWaRD')
			legendStr(M+1) = 'Schiske';
			legend(legendStr);
			xlabel('BSNR (dB)')
			ylabel('Relative l^2 error')
		
	otherwise
		
end
xlim([min(levels) max(levels)])
