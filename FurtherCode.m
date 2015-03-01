H0 = totalAllP;

log_liklyhood_bound = [];
log_liklyhood_rand = [];
log_liklyhood_all = [];
for CH = 1:22

log_liklyhood = hists{CH}' * log(H0);
boundBlocksTmp = bounds{CH};
if boundBlocksTmp(end,2) > size(hists{CH},2)
	boundBlocksTmp(boundBlocksTmp/40000 > size(hists{CH},2)) = 0;
	boundBlocksTmp = boundBlocksTmp(all(boundBlocksTmp,2),:);
end
boundBlocks = union(boundBlocksTmp(:,1),boundBlocksTmp(:,2))/40000;
allBlocks = (1:size(log_liklyhood,1))';
randBlocks =  randperm(size(allBlocks,1),size(boundBlocks,1))';

log_liklyhood_bound = cat(1,log_liklyhood_bound , log_liklyhood(boundBlocks));
log_liklyhood_rand = cat(1,log_liklyhood_rand , log_liklyhood(randBlocks));
log_liklyhood_all = cat(1,log_liklyhood_all , log_liklyhood(allBlocks));
end

log_liklyhood_bound = log_liklyhood_bound(log_liklyhood_bound<-2.5*1e5);
log_liklyhood_rand = log_liklyhood_rand(log_liklyhood_rand<-2.5*1e5);
log_liklyhood_all = log_liklyhood_all(log_liklyhood_all<-2.5*1e5);

[log_liklyhood_bound_hist, bound_edges] = histcounts(log_liklyhood_bound,40);
[log_liklyhood_rand_hist, rand_edges] = histcounts(log_liklyhood_rand,bound_edges);
[log_liklyhood_all_hist, all_edges] = histcounts(log_liklyhood_all,bound_edges);

log_liklyhood_bound_hist = log_liklyhood_bound_hist./sum(log_liklyhood_bound_hist);
log_liklyhood_rand_hist = log_liklyhood_rand_hist./sum(log_liklyhood_rand_hist);
log_liklyhood_all_hist = log_liklyhood_all_hist./sum(log_liklyhood_all_hist);

bound_x = (bound_edges(1:end-1) + bound_edges(2:end))/2;
rand_x = (rand_edges(1:end-1) + rand_edges(2:end))/2;
all_x = (all_edges(1:end-1) + all_edges(2:end))/2;

plot(bound_x,log_liklyhood_bound_hist, rand_x,log_liklyhood_rand_hist, all_x,log_liklyhood_all_hist);
legend('Bound','Rand','All','Location','northwest');
%figure;plot(bound_x,log_liklyhood_bound_hist);
%figure;plot(rand_x,log_liklyhood_rand_hist);
%figure;plot(all_x,log_liklyhood_all_hist);
