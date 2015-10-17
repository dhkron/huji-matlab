function antidg = antidiag(mat, dpos)
	% DPOS = location along diagonal
	% This antidig goes through DPOS
	% Offset = 0 is normal
	% Offset = a_size for invisible antidiag
	a_size = size(mat,1);
	start = (a_size+1)*dpos - a_size;
	step = a_size-1;
	if dpos == a_size/2
		elemcount = a_size/2;
	elseif dpos > a_size/2
		elemcount = a_size-dpos+1;	
	else
		elemcount = dpos;	
	end
	antidg = mat(start:step:end);
	antidg = antidg(1:elemcount);
	return
end
