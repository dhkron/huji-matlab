%Compare avg of all blocks to:
% 1) boundery blocks
% 2) same amount of random blocks
%Then, compare the bound results for different chr's or subsets of bound.

function [diffRandom, diffBound] = MotivComp(blockHistVector, bounds, blockSize)

%Infer variables
numBlocks = size(blockHistVector,2);
words = 1:size(blockHistVector,1);

%Calculate block numbers
allBlocks = 1:numBlocks;
boundBlocksStart = bounds(:,1)'/40000 + 1; %offset 40000 = block#2
boundBlocksEnd = bounds(:,2)'/40000 + 1;
randBlocks = randperm(numBlocks,size(boundBlocksStart,2)); %Select random collection of blocks, with count euqal to bounds'

yrng = [-5e-4 5e-4];

allHist = MotivHist(allBlocks, blockVector, wordLen, wordBase);
randHist = MotivHist(randBlocks, blockVector, wordLen, wordBase);
boundHist = MotivHist(boundBlocks, blockVector, wordLen, wordBase);

%figure('Name','rand');plot(1:1024, randHist-allHist); ylim(yrng);

%accuracy = 1-[0.1 0.25 0.5 0.75 0.9 1];
%figure;
%hold on;
%for acc = accuracy
%	boundBlocksR = boundBlocks(rand(1,size(boundBlocks,2))>=acc);
%	boundHist = MotivHist(boundBlocksR, blockVector, wordLen, wordBase);
%	figure('Name',sprintf('hist %.2f',acc));
%	plot(1:1024, boundHist-allHist); ylim(yrng);
%end
%legend(arrayfun( @(x) sprintf('hist %.2f',x), accuracy, 'Unif',false));
%hold off;

diffRandom = randHist-allHist;
diffBound = boundHist-allHist;

plot(1:wordBase^wordLen, diffBound./diffRandom);

figure
subplot(2,1,1);
plot(1:wordBase^wordLen, diffRandom);
ylim(yrng);
legend('rand');

subplot(2,1,2);
plot(1:wordBase^wordLen, diffBound);
ylim(yrng);
legend('bound');

end
