function [a,a_norm] = RegularLoader(chrNum)
	dataPath = sprintf('/cs/cbio/gil/raw/Liz/WT_chr%d.txt',chrNum);
	normPath = sprintf('HiC-CSV-Matrices/normalization_vector_LIZ_WT_chr%d.mat', chrNum);	
	%dataPath = sprintf('/cs/cbio/gil/raw/Liz/KO_chr%d.txt',chrNum);
	%normPath = sprintf('HiC-CSV-Matrices/normalization_vector_LIZ_chr%d.mat', chrNum);	

	a = load(dataPath);

	if exist(normPath,'file')
		%fprintf('Norm vector found');
		load(normPath,'b');
	else
		fprintf('Recalculating norm vector');
		b = Normalize(a);
		fprintf('Saving to dist');
		save(normPath,'b');
	end

	a_norm = a ./ (b'*b);
	
	a_clean = a_norm;
	a_clean(isnan(a_clean)) = 0;
	a_sums = sum(a_clean);
	med = median(a_sums(a_sums>0));
	a_norm = a_norm/med;

	%Normalize
	%% The problematic cells are killed by b, having inf in their location!
	%% Any following code MUST drop zero value cells.
	%% This also kills other zeros cells
	% ( sum(a)==0) ) == (b==Inf)
	% These are the indices of interest
	% ind = (b!=Inf)
	% a_clean = a(ind,ind)
	% a_orig(ind,ind) = a_clean
	% When removing, TAD structure is curios as well...
end
