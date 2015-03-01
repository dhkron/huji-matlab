r1 = 1;
r2 = 300;

rleft = ceil(r1/sqrt(2));
rright = floor(r2/sqrt(2));

dim = rright-rleft+1;

disp(sprintf('range %d %d dim %d',rleft,rright,dim))
disp(sprintf('invalid blocks %d',sum(sum(blockVector(:,rleft:rright)>=0)<40000)))

a = zeros(dim,dim);
for i = 1:dim
	for j = 1:dim
	[~,vec1] = GenomeHist(i+rleft-1,blockVector,5,4);
	[~,vec2] = GenomeHist(j+rleft-1,blockVector,5,4);
	a(i,j) = corr(vec1',vec2');
	end
end
figure;
imagesc(a);
colorbar;
axis equal;
