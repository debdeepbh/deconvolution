function [out] = demo_sig_perf() 

nums = 1:6;

levels = 1:30;

scaling = 1;

schiske_err = zeros(size(nums), size(levels));
schiskeforwd_err = zeros(size(nums), size(levels));

for num=nums
	ind = 1;
	for level=levels	
		[schiske_err(num, ind), schiskeforwd_err(num, ind)] = demo_schiske_vs_schiskeforwd(num, level, scaling);
		ind = ind+1;
	end
	%figure(testvec_num)
	%plot(levels, schiske_err, levels, schiskeforwd_err);
	%legend('Schiske', 'Multichannel ForWaRD');
	%xlim([min(levels) max(levels)])
end

for num = nums
 	plot(levels, schiskeforwd_err(num,:));
	hold on
end









