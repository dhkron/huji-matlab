%This one defines a function that draw the motiv hist for given blocks

% wordLen is usually 5
% wordBase is usually 4
function [aHist] = MotivHist(blockNumbers, blockVector, wordLen, wordBase)
	aHist = zeros(1,wordBase^wordLen);
	count = 0;
	for blockNum = blockNumbers
		[~,tmphist] = GenomeHist(blockNum,blockVector,wordLen,wordBase);
		if sum(tmphist) ~= 39995
			count = count + 1;
		end
		aHist = aHist + tmphist;
	end
	count,
	aHist = aHist / sum(aHist);%size(blockNumbers,2);
end
