function LlrHeatmap(a,i,j,h)
	%Returns 2D matirx L
	% L(i,j) = sum(llr(heatmap(range i to j, antidiagonally)))
	llr_sum = 0;
	b = a(1:end-1,2:end);
	for pos = i:j
		antidg = Antidiag(a,pos);
		antidg = antidg(1:min(end,h))
		llr_sum = llr_sum + sum(antidg);
		antidg = Antidiag(b,pos);
		antidg = antidg(1:min(end,h))
		llr_sum = llr_sum + sum(antidg);
	end
end
