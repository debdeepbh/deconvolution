% get the filename of the impulse response of the antenna, given the antenna absolute index
function filename = getantenna(num)
% num from 0 to 47

polarity = 'H';

T = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
M = [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
B = [32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47];

if (length(find(T == num)) == 1)
	location = 'T';
	ind = find(T==num);
elseif (length(find(M == num)) == 1)
	location = 'M';
	ind = find(M==num);
elseif (length(find(B == num)) == 1)
	location = 'B';
	ind = find(B==num);
end

filename = strcat(sprintf('%02d',ind),location,polarity,'.imp');



