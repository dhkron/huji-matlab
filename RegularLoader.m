function [a,a_norm] = RegularLoader(chrNum)
	dataPath = sprintf('HiC-CSV-Matrices/IMR90_chr%d_40k',chrNum);
	normPath = sprintf('HiC-CSV-Matrices/normalization_vector_IMR90_chr%d_40k.mat', chrNum);	

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
	med = median(a_sums(a_sums>0)); %Median is equal to everyone here
	a_norm = a_norm/med; %Normalize row sums to 1
	%a_norm = a_norm/sqrt(sum(sum(b'*b~=inf))/numel(b)^2);

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
