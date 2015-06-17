function [a_llr,a_pdt,a_pdb] = PosteriorHeatmap(a,MIN_DIAG,MAX_DIAG,probBg,meansBg,sigmaBg,probIn,meansIn,sigmaIn,box)

a_size = size(a,1);
a_pdt = zeros(a_size);
a_pdb = zeros(a_size);
for i = MIN_DIAG:MAX_DIAG
	start_pos = a_size*i + 1;
	diag_elem = GetDiag(a,i,0,0,1,1); %Still maps zero to zero

	pdt_diag = normpdf(diag_elem,meansIn(i),sigmaIn(i))*probIn(i);
	pdt_diag(diag_elem==0) = 0;
	a_pdt(start_pos:a_size+1:end) = pdt_diag;

	pdb_diag = normpdf(diag_elem,meansBg(i),sigmaBg(i))*probBg(i);
	pdb_diag(diag_elem==0) = 0;
	a_pdb(start_pos:a_size+1:end) = pdb_diag;
end

	a_pdsum = a_pdt + a_pdb;
	a_pdt = a_pdt ./ a_pdsum;
	a_pdb = a_pdb ./ a_pdsum;
	a_pdt(isnan(a_pdt)) = 0;
	a_pdb(isnan(a_pdt)) = 0;
	a_llr = log2(a_pdt)-log2(a_pdb); %rename LLR
	a_llr(isnan(a_llr)) = 0;

	if exist('box','var')
		a_pdt = a_pdt(box,box);
		a_pdb = a_pdb(box,box);
		a_llr = a_llr(box,box);
	end
end
