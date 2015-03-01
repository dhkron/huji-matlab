%Define constants
function [blockVector, NGrams, expBase, cmap] = Genome(chrNumber)
dataFile = sprintf('chr%d',chrNumber);
blockSize = 40000;
charsPerRow = 50;
NGrams = 5; % 5, wordlen
expBase = 4; %a,c,t,g
cmap = zeros(1,116)-expBase^(NGrams+1);
%cmap('nN') = -expBase^(NGrams+1);
cmap('aA') = 0;
cmap('cC') = 1;
cmap('gG') = 2;
cmap('tT') = 3;

disp(sprintf('Data %s Ngrams %d Base %d',dataFile,NGrams,expBase))

%Open files
fid = fopen(sprintf('/cs/cbio/gil/bowtie/%s.fa',dataFile),'r');
C = textscan(fid, '%s', 'delimiter', '\n','headerlines',1);

%Generate genome-wide data vector
%After this code, blockVector contains blocks of DNA as numbers
lastLine = size(C{1},1);
lastLineSize = size(C{1}{lastLine},2);
blockVector = zeros(blockSize, ceil(size(C{1},1)*charsPerRow/blockSize));
for lineNum = 1:size(C{1},1)
	currentBlock = floor((lineNum-1)*charsPerRow/blockSize);
	currentPos = mod((lineNum-1)*charsPerRow,blockSize);
	realBlock = currentBlock + 1;
	realPos = currentPos + 1;
	if lineNum == lastLine
		blockVector(realPos:(realPos+lastLineSize-1),realBlock) = cmap(C{1}{lineNum});
	else
		blockVector(realPos:(realPos+charsPerRow-1),realBlock) = cmap(C{1}{lineNum});
	end
end
end
