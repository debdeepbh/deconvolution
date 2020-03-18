function [testvec] = testvec_gen(n) 
testvec = zeros(1,512);
switch n
	case 1	% Two chirps
		for j = 129:256
			n = j-1;
			testvec(j) = sin (((abs(n - 128))^1.7) / 128);
		end
		for j = 385:448
			n = j-1;
			testvec(j) = sin (((abs(n - 128))^2) / 128);
		end

	case 2	% boxes
		for j = 32:95
			testvec(j+1) =	1;
		end
		for j = 132:259
			testvec(j+1) =	2;
		end
		for j = 325:415
			testvec(j+1) =	4;
		end
	case 3	% sawtooth/triangle
		for j = 0:63
			testvec(j+1) = 1 - j/64;
		end
		for j = 256:319
			testvec(j+1) = 5 - j/64;
		end
		for j = 64:255
			testvec(j+1) =0;
		end
		for j = 320:511
			testvec(j+1) =0;
		end
	case 4	% chirp
		for j = 0:511
			testvec(j+1) = sin(j^(1.5)/64);
		end

	case 5	% smooth impulsive
		for j = 0:511 
		       testvec(j+1) = (j - 256)*exp(-(j-256)^2/512);
		end
	case 6	% narrow impulsive
		for j = 0:511 
		       testvec(j+1) = (j - 256)*exp(-(j-256)^2/32);
		end
	otherwise
		error('Wrong argument: %d', n)
		
end

