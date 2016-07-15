function BreakMatrix(a,psize,overlap,pathbase)
	v = [1 psize]
	sz = size(a,2);
	while v(end) < sz
		v(end+1,:) = v(end,:)+psize-overlap;
	end
	v(end,2)=sz;
	v(end,1)=sz-psize+1;
	for i = v'
		box = i(1):i(2);
		M = a(box,box);
		path = sprintf('%s.%d_%d.txt',pathbase,i(1),i(2));
		Log(['Saving ' path])
		dlmwrite(path,M,'precision',20,'delimiter','\t')
		Log()
	end
end
