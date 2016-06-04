function [a,a_norm] = LibermanLoader(chrNum,res)
	dirOfMight = '/cs/cbio/gil/raw/GSE63525_5K_Liberman/GM12878_combined/%sb_resolution_intrachromosomal/chr%d';
	mapq1 = 'MAPQG0'; %More data here, less accurate
	mapq2 = 'MAPQGE30';
	filePrefix = 'chr%d_%sb';
	fileTypeData = 'RAWobserved';
	fileTypeNorm = 'KRnorm';

	resKb = res/1000;
	if resKb < 1000
		resPrnt = sprintf('%dk',resKb);
	else 
		resPrnt = sprintf('%dm',resKb/1000);
	end

	realDir = sprintf(dirOfMight,resPrnt,chrNum); %res = 10
	realFilePrefix = sprintf(filePrefix,chrNum,resPrnt);
	dataPath = sprintf('%s/%s/%s.%s',realDir,mapq2,realFilePrefix,fileTypeData);
	normPath = sprintf('%s/%s/%s.%s',realDir,mapq2,realFilePrefix,fileTypeNorm);

	a_raw = load(dataPath);
	a_raw(:,1) = a_raw(:,1)/res+1; %So that first res-base block is 1
	a_raw(:,2) = a_raw(:,2)/res+1; %So that first res-base block is 1
	a = spconvert(a_raw);
	
	b = load(normPath);

	sza = size(a);
	szn = size(b,1);
	assert(szn > 1);
	if sza(1) == sza(2) && sza(1) == szn
		%Nothing to do here, square matrix at perfect size
	else
		a(szn,szn) = 0;
	end
	a_norm = a ./ (b*b');
	
	%Find median
	%a_clean = a_norm;
	%a_clean(isnan(a_clean)) = 0;
	%a_sums = sum(a_clean);
	%med = median(a_sums(a_sums>0));
	%a_norm = a_norm/med;
end
