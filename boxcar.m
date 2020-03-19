function [out] = boxcar(supp, N) 
% outputs a characteristic function with support size supp
out = zeros(1, N);
out(1:supp) = 1;

