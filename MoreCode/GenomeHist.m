%Outputs an histogram according to blocknumber
function [wordBlock, histVec] = GenomeHist(targetBlockNum,blockVector,wordlen,wordbase)
	wordBlock = GenomeToWordVec(blockVector(:,targetBlockNum), wordlen, wordbase);
	histVec = histc(wordBlock,0:(wordbase^(wordlen)-1))';
end

%expCurrent = 1;
%targetBlock = blockVector(:,targetBlockNum);
%wordBlock = zeros(size(targetBlock,1)-NGrams+1,1); 
%for i = NGrams:-1:1
%	wordBlock = wordBlock + expCurrent*targetBlock(i:end-NGrams+i);
%	expCurrent = expCurrent*expBase;
%end
