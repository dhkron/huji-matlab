%Compare avg of all blocks to:
% 1) boundery blocks
% 2) same amount of random blocks
%Then, compare the bound results for different chr's or subsets of bound.

function [bData] = MotivComp(blockHistVector, bounds, blockSize, plotOne, plotAll)

%Infer&define variables
numBlocks = size(blockHistVector,2);
words = 0:size(blockHistVector,1)-1;
yrng = [-5e-4 5e-4];

%
if bounds(end,2) > size(blockHistVector,2) 
	disp 'Bounds too large, truncating'
	bounds(bounds/blockSize > size(blockHistVector,2)) = 0;
	bounds = bounds(all(bounds,2),:);
end

%Calculate block numbers
allBlocks = 1:numBlocks;
boundStart = bounds(:,1)'/blockSize + 1; %offset 40000 = block#2
boundEnd = bounds(:,2)'/blockSize + 1; %offset 40000 = block#2
boundAll = union(boundStart,boundEnd);
boundStartUniq = setdiff(boundStart,boundEnd);
boundEndUniq = setdiff(boundEnd,boundStart);
boundCommon = intersect(boundStart,boundEnd);
randBlocks = randperm(numBlocks,size(boundAll,2)); %Select random collection of blocks, with count equal to bounds'

%Arrange block numbers in data structure
% 1	2	 3	4	 5
% name	blockset hist 	normhist (normhist-allhist)
bData = {};
bData{1} = {'all',allBlocks};
bData{2} = {'rand',randBlocks};
bData{3} = {'blockAll',boundAll};
bData{4} = {'blockCommon',boundCommon};

%bData{3} = {'blockStartUniq',boundStartUniq};
%bData{4} = {'blockEndUniq',boundEndUniq};
%bData{7} = {'blockStart',boundStart};
%bData{8} = {'blockEnd',boundEnd};

%bData{5} = {'blockend3',boundBlocksEnd+3};
%bData{6} = {'blockend9',boundBlocksEnd+9};
%bData{7} = {'blockend15',boundBlocksEnd+15};

%Calculate normalized hist for each block set
for i = 1:size(bData,2)
	blockSet = bData{i}{2};
	sumHist = sum(blockHistVector(:,blockSet),2);
	normHist = sumHist/sum(sumHist);
	bData{i}{3} = sumHist;
	bData{i}{4} = normHist;
	bData{i}{5} = normHist-bData{1}{4};
end

%Plot all hists
if exist('plotAll','var') && plotAll
	for i = 2:size(bData,2)
		figure('name',bData{i}{1});plot(words, bData{i}{4})
	end
end

%Chi
% N * null hypo
E = sum(bData{2}{3}) * bData{1}{4};
Chi = sum(((bData{2}{3}-E).^2)./E);
1-chi2cdf(Chi,1024);

E = sum(bData{3}{3}) * bData{1}{4};
Chi= sum(((bData{3}{3}-E).^2)./E);
1-chi2cdf(Chi,1024);

E = sum(bData{4}{3}) * bData{1}{4};
Chi= sum(((bData{4}{3}-E).^2)./E);
1-chi2cdf(Chi,1024);

%plot(words, diffBound./diffRandom);

if exist('plotOne','var') && plotOne
	figure('name','block comparison')
	numPlots = size(bData,2);

	for i=2:numPlots
		subplot(numPlots,1,i-1);
		plot(words, bData{i}{5});
		ylim(yrng);
		legend(bData{i}{1});
	end
end

end
