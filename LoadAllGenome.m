wordlen = 8;
if ~exist('hists') || size(hists,2)~=22
	hists = {};
	for i=1:22
		hists{i} = LoadGenomeToHists(i,40000,wordlen,4);
	end
end
if ~exist('bounds') || size(bounds,2)~=22
	bounds = {};
	for i=1:22
		bounds{i} = load(sprintf('TADs/bounds.IMR90.chr%d',i));
	end
end
if ~exist('bDatas') || size(bDatas,2)~=22
	bDatas = {};
	for i=1:22
		bDatas{i} = MotivComp(hists{i},bounds{i},40000);
	end
end

indexRawHist = 3;
totalAllHist = bDatas{1}{1}{indexRawHist};
totalRandHist = bDatas{1}{2}{indexRawHist};
totalBoundHist = bDatas{1}{3}{indexRawHist};
for i=1:22
	totalAllHist = totalAllHist + bDatas{i}{1}{indexRawHist};
	totalRandHist = totalRandHist + bDatas{i}{2}{indexRawHist};
	totalBoundHist = totalBoundHist + bDatas{i}{3}{indexRawHist};
end

totalAllP = totalAllHist/sum(totalAllHist);
totalRandP = totalRandHist/sum(totalRandHist);
totalBoundP = totalBoundHist/sum(totalBoundHist);

