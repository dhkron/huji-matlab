function [a_lpr,a_pdt,a_pdb] = LiklihoodHeatmap(a,MIN_DIAG,MAX_DIAG,meansBg,sigmaBg,meansIn,sigmaIn)

probBg = ones(size(meansBg));
probIn = ones(size(meansIn));
[a_lpr,a_pdt,a_pdb] = PosteriorHeatmap(a,MIN_DIAG,MAX_DIAG,probBg,meansBg,sigmaBg,probIn,meansIn,sigmaIn,true);

end
