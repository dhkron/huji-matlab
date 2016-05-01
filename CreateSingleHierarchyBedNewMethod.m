function [a,a_llr,s] = CreateSingleHierarchyBedNewMethod(fMat,fNewDomains,fBgModel,prefix,res,chr,box,bedPath,figPath)
	hasBox = ( exist('box','var') && numel(box)>1 );
	hasBed = ( exist('bedPath','var') && numel(bedPath)>0 );
	hasFig = ( exist('figPath','var') && numel(figPath)>0 );

	if ischar(res)
		res = str2num(res);
	end;
	if ischar(chr)
		chr = str2num(chr);
	end;

	Log('Loading files');
	a = load(fMat);
	bnd = load(fNewDomains);
	boundaries = (bnd(2:end,1)'/res)+1;
	M = load(fBgModel);
	bgModel = M.meansBgDxn;
	Log();

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

	Log('Constructing hierarchy');
	RegressMerger(boundaries,a,box(1),box(end),res,bgModel,chr,bedPath);
	%TadTree(d1,s,find(d3),box(1),box(end),a_llr,-a_llr,res,chr,bedPath);
	Log();

	view(-45,90);

	title(sprintf('Hierarchy TAD numbered nesting - %s chr%d [blocks %d-%d]',prefix,chr,box(1),box(end)));

	g=gcf;
	if hasFig
		SaveFigure(g,figPath);
		%path = sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/Hierarchies_%s_Chr%d_Full.png',prefix,chr);
	end
	close(gcf);
end

