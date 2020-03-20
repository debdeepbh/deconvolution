function [opt_sc, ratio_above ] = getopt_j(B_j, L, r, noiseax, w_ori_j, rootmethod) 
% computing the optimum scaling paramer
% for which ratio is equal to the parameter
% alpha_j = getratio_above_noise(alpha_j, B_j, L, r, noiseax, w_ori_j) 

switch rootmethod
case 'bisec'
	% initial guess
	a = 0.001;
	%a = 0.02;
	b = 0.999;

	% error
	errorl = 0.01;

	while (abs(a-b)>errorl)
		%alpha(i) = a;
		opt_sc = a;
		printf('.')

		%[w, ratiounthr, thrvec] = wienforwd(z,k,type,p,sigma,alpha, rho,method);
		%fa = ratiounthr(i);
		fa = getratio_above_noise_schiske(opt_sc, B_j, L, r, noiseax, w_ori_j);

		%alpha(i) = b;
		opt_sc = b;
		%[w, ratiounthr, thrvec] = wienforwd(z,k,type,p,sigma,alpha, rho,method);
		%fb = ratiounthr(i);
		fb = getratio_above_noise_schiske(opt_sc, B_j, L, r, noiseax, w_ori_j);


		if ((a-fa)>=0 && (b - fb)<=0)
			c = (a+b)/2;
			%alpha(i) = c;
			opt_sc = c;

			%[w, ratiounthr, thrvec] = wienforwd(z,k,type,p,sigma,alpha, rho,method);
			%fc = ratiounthr(i);
			fc = getratio_above_noise_schiske(opt_sc, B_j, L, r, noiseax, w_ori_j);

			if ((c -fc) >0)
				a = c;
			else
				b = c;
			end
		elseif ((a-fa)<=0 && (b - fb)>=0)
			c = (a+b)/2;
			%alpha(i) = c;
			opt_sc = c;
			%[w, ratiounthr, thrvec] = wienforwd(z,k,type,p,sigma,alpha, rho,method);
			%fc = ratiounthr(i);
			fc = getratio_above_noise_schiske(opt_sc, B_j, L, r, noiseax, w_ori_j);

			if ((c -fc) >0)
				b = c;
			else
				a = c;
			end
		else
			c = 1;
			fc = getratio_above_noise_schiske(1, B_j, L, r, noiseax, w_ori_j);
			disp('no solution in between');
			break;
		end
	end
	opt_sc = c;
	ratio_above = fc;

case 'search'
	sc_initial = 0.01;
	sc_step    = 0.001;
	sc_end     = 1;

	% brute force
	min_val = 10e4;
	opt_sc = 1;
	ratio_above = 1;

	for sc_val=sc_initial:sc_step:sc_end
		ratio_above = getratio_above_noise_schiske(sc_val, B_j, L, r, noiseax, w_ori_j);
		% compute the difference
		y = abs(ratio_above  - sc_val);

		% compute the running minimum
		min_val = min(min_val, y);

		% if running minimum is the current val, store the scaling value
		if min_val == y
			opt_sc = sc_val;
			opt_ratio = ratio_above;
		end
	end

	%% If the minimizer does not exist in [0,1]
	%% we prefer to choose opt_sc =1 to be close to Schiske
	if opt_sc == sc_step
		opt_sc = 1;
		opt_ratio = getratio_above_noise_schiske(1, B_j, L, r, noiseax, w_ori_j);
	end       

case 'noreg'

	% original signal coeffs above the standard noise level (with scaling 1)
	opt_sc = getratio_above_noise_schiske(1, B_j, L, r, noiseax, w_ori_j);
	ratio_above = opt_sc;

otherwise
	disp('wrong rootfinding method')
end






