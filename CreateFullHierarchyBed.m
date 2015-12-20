function [a,a_llr,s] = CreateFullHierarchyBed(prefix,res)
for chr = 1:19
	fprintf('Starting chr %d\r\n',chr);
	template = '/cs/cbio/gil/htad-chain/output/%s.%s.chr%d.txt';
	a = load(sprintf(template,prefix,'matrix',chr));
	a_llr = load(sprintf(template,prefix,'llr',chr));
	s = load(sprintf(template,prefix,'supersum',chr));
	box = 1:size(a,2);
	DisplayHeatmap(log2(a+1),[0,6],box,'red');[d1 d2 d3]=DynProgTAD(s,box(1),box(end),0);TadTree(d1,s,find(d3),box(1),box(end),a_llr,-a_llr,res,chr,prefix);

	title(sprintf('Hierarchy TAD numbered nesting - %s chr%d [blocks %d-%d]',prefix,chr,box(1),box(end)));

	g=gcf;
	path = sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/Hierarchies_%s_Chr%d_Full.png',prefix,chr);
	SaveFigure(g,path);
	close(g);
end

