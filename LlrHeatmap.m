function llrmap = LlrHeatmap(a,h)
	%Returns 2D matirx L
	% L(i,j) = sum(llr(heatmap(range i to j, antidiagonally)))
	llrmap = zeros(size(a));
	llrstt = zeros(size(a));
	a_size = size(a,1);
	b = a(1:end-1,2:end);
	for d=0:a_size %Diagonal number
		for p=1:a_size-d
			x = p;
			y = d + p;
			llrmap(x,y) = LlrHeatmapHelper(a,b,x,y,h);
		end
	end

function psum = LlrHeatmapHelper(a,b,i,j,h)
	%Rewrite this whole thing. Very unefficient.
	%Built it bottom up WITHOUT ifs.
	if i==j
		antidg = Antidiag(a,i);
		antidg = antidg(1:min(end,h));
		psum = sum(antidg);
	elseif i<j
		mid = floor((i+j)/2);
		psum1 = llrmap(i,mid);
		psum2 = llrmap(mid+1,j);
		antidg = Antidiag(b,mid);
		antidg = antidg(1:min(end,h-1));
		psum = psum1+psum2+sum(antidg);
	end
end

end
