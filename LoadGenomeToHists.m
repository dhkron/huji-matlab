%This function is used to load a .fa file of a chromosome,
% and store it as a vector of histogram blocks.
%Stages:
% 1. Load genome file
% 2. Convert nt to numbers
% 3. Calculate words (
% 4. Reshape to blocks (40k typical)
% 5. For each block, calculate histogram (of positive words only)
% 6. Store histograms
function [blockHistVector, cmap] = LoadGenomeToHists(chrNumber, blockSize, wordlen, wordbase)
%Constants
dataFile = sprintf('chr%s',chrNumber);
cmap = zeros(1,116)-wordbase^(wordlen+1); %Ensure negativity of words containing bad chars
cmap('aA') = 0;
cmap('cC') = 1;
cmap('gG') = 2;
cmap('tT') = 3;

disp(sprintf('Data %s Ngrams %d Base %d Blocksize %d',dataFile,wordlen, wordbase, blockSize));

totalTime = tic;
%Open & read file
fprintf('Open & read... ');tic;
fid = fopen(sprintf('/cs/cbio/gil/bowtie/%s.fa',dataFile),'r');
C = textscan(fid, '%s', 'delimiter', '\n','headerlines',1);
fprintf('Done in %3.3f seconds\r\n',toc);
%Merge the string
fprintf('Merge string... ');tic;
longString = sprintf('%s',C{1}{:});
fprintf('Done in %3.3f seconds\r\n',toc);
%Convert to numbers
fprintf('Convert to numbers... ');tic;
longVector = cmap(longString);%Row vector
fprintf('Done in %3.3f seconds\r\n',toc);
%Convert to words (as ints)
fprintf('Convert to words... ');tic;
longWordVector = GenomeToWordVec(longVector, wordlen, wordbase);%Row vector
fprintf('Done in %3.3f seconds\r\n',toc);
%Trim and convert to blocks
fprintf('Reshaping... ');tic;
blockVector = reshape(longWordVector(1:end-mod(end,blockSize)),blockSize,[]);
fprintf('Done in %3.3f seconds\r\n',toc);
%Get histograms
fprintf('Calculating histograms... ');tic;
blockHistVector = histc(blockVector,0:(wordbase^(wordlen)-1));
fprintf('Done in %3.3f seconds\r\n',toc);
%Total time
fprintf('Finished in %3.3f seconds\r\n',toc(totalTime));

end

%Some integ testing snippets
% blocknum = 30
% sprintf('%d',cmap(C{1}{(blocknum-1)*40000/50+1}(2:6)))
% dec2base(x(2,blocknum),4,5)
