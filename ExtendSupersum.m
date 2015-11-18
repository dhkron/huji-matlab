%
% Add to s the value at x,y - and also any things one the way
% Assum that x < y (top pyramid)
%
function v = ExtendSupersum(x,y,a_t,a_b)
	% L(i,j) = sum( llr(pyramid is tad) + llr(sky are background) )

	%The sum of the TAD pyramid at that location.
	%Since a_t is zero in the lower diag, should be ok to take the box
	pyrIsTad = sum(sum( a_t(x:y,x:y) ));

	%The sum of background of the sky area
	%It is the sum of the square minus the sum of the columns
	%More easy to fliplr the matrix, because then diag returns the right result.
	skyAreBlue = -sum(sum( a_b(x:y,x:y) ));

	a_bg = fliplr(a_b);
	a_size = size(a_t,1);
	dg_s = a_size+1-2*y;
	dg_e = a_size+1-2*x;
	for i=dg_s:dg_e
		skyAreBlue = skyAreBlue + sum(diag(a_bg,i));
	end

	v = pyrIsTad + skyAreBlue;
end

