test_count = 1000
test_size = 1e6
disp 'Testing translation to words'
gene1=ceil(rand(1,test_size)*4);%Use 1-4 indexing
disp 'sub2ind'
tic;
for i=1:test_count
	out1 = sub2ind([4 4 4 4 4],gene1(5:end),gene1(4:end-1),gene1(3:end-2),gene1(2:end-3),gene1(1:end-4));
end
toc;

gene2=gene1-1;%Using 0-3 indexing
disp 'gil func'
tic;
for i=1:test_count
	out2 = zeros(1,size(gene2,2)-4);
	expc = 1;
	for k=5:-1:1
		out2 = out2 + expc * gene2(k:end-5+k);
		expc = expc * 4;
	end
end
toc;

disp 'Should be zero: ',sum(out1-out2-1)
