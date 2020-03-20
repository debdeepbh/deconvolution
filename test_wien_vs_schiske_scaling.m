function [wien_err, schiske_err] = test_wien_vs_schiske_scaling(testvec_num, noise_level) 

sclist = 0.1:0.1:1.5;
i=1;
for sc=sclist;
	[wien_err(i), schiske_err(i)] = demo_wien_vs_schiske_scaling(testvec_num, noise_level, sc);
	i = i + 1;
end

plot(sclist, wien_err, sclist, schiske_err, sclist, min(wien_err) * (zeros(size(sclist)) + 1), sclist, min(schiske_err) * (zeros(size(sclist)) + 1));
legend('wien err', 'schiske err')
xlim([min(sclist) max(sclist)]);


	

