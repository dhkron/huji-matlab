function [a,a_llr,s] = CreateSingleHierarchyBed(fMat,fLLR,fSupersum,prefix,res,chr,box,bedPath,figPath)
	hasBox = ( exist('box','var') && numel(box)>1 );
	hasBed = ( exist('bedPath','var') && numel(bedPath)>0 );
	hasFig = ( exist('figPath','var') && numel(figPath)>0 );

	a = load(fMat);
	a_llr = load(fLLR);
	s = load(fSupersum);

	if ischar(res)
		res = str2num(res);
	end;
	if ischar(chr)
		chr = str2num(chr);
	end;

	if ~hasBox
		box = 1:size(a,2);
	end

	if ~hasBed
		bedPath = '';
	end
	if ~hasFig
		figPath = '';
	end

	DisplayHeatmap(log2(a+1),[0,6],box,'red');
	[d1 d2 d3]=DynProgTAD(s,box(1),box(end),0);
	TadTree(d1,s,find(d3),box(1),box(end),a_llr,-a_llr,res,chr,bedPath);

	view(-45,90);

	title(sprintf('Hierarchy TAD numbered nesting - %s chr%d [blocks %d-%d]',prefix,chr,box(1),box(end)));

	g=gcf;
	if hasFig
		SaveFigure(g,figPath);
		%path = sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/Hierarchies_%s_Chr%d_Full.png',prefix,chr);
	end
	close(gcf);
end

