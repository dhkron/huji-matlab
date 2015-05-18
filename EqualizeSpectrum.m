function a_eq = EqualizeSpectrum(a, MIN_DIAG, MAX_DIAG, plotStuff)
	a_size = size(a,1);
	a_eq = zeros(a_size);
 	refdiag = diag(a,MINDIAG);
	refsort = sort(refdiag);
	for i = MIN_DIAG+1:MAX_DIAG
	 	b = diag(a,i);
		bsorted = sort(b);
		f = i-MIN_DIAG-1;
		%c1 = corr(refdiag(1:end-f),b);
		%c2 = corr(refsort(1:end-f),bsorted);

		% Estimate linearly
		p = polyfit(bsorted, refsort(1:end-f),1);
		polydata = polyval(p,bsorted);
		% Fill data in a_eq
		start_pos = a_size*i + 1;
		a_eq(start_pos:a_size+1:end) = polydata;
	end
	if exist('plotStuff','var') && plotStuff
		figure;
		for i = MIN_DIAG:MAX_DIAG
			dgn = diag(a,i);
			dgn = sort(dgn);
			plot(dgn)
			hold on;
		end
		title('Spectrum of diagonals');
		hold off;
	end
end
