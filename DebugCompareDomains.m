function [] = DebugCompareDomains(fmat,fsum,fdom,res,box)
	a = load(fmat);
	s = load(fsum);
	d = load(fdom);

	if ischar(res)
		res = str2num(res);
	end

	if exist('box','var') && numel(box)>1
		%All good, box is here
	else
		box = 1:size(a,2);
	end

	cost = 0;

	%My domains
	DisplayHeatmap(log(a),[0 10],box,'red');
	[d1 d2 d3] = DynProgTAD(s,box(1),box(end),cost);
	view(-45,90);
	title('My Method');

	%Dixon
	DebugDixonDomains(a,d,res,box);
	caxis([0 10]);
	view(-45,90);
	title('Dixon Method w 5 iter');
end
