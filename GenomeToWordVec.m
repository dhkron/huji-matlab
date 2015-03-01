%This function converts a gene vector into a word vector.
% - It assumes gene vector is composed of numbers 0..wordbase-1
% - A user should make sure that any numbers outside this range
% will be very negative. That is, a string of the form 'tttttttttn' will 
% always be negative. Thus, n should be mapped to -minus wordbase^(wordlen+1)
% - 'AAAAT' is translated to 00003, keeping the normal order in mind
% - Ah, and it's all row vectors. Might cause problems in existing scripts
function [wordVector] = GenomeToWordVec(geneVector,wordlen,wordbase)
	expCurrent = 1;
	wordVector = zeros(1,size(geneVector,2)-wordlen+1); 
	for i = wordlen:-1:1
		wordVector = wordVector+expCurrent*geneVector(i:end-wordlen+i);
		expCurrent = expCurrent*wordbase;
	end
end

