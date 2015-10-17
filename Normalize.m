%
% Normalizes Hi-C data
%
% Assumption: Data is more or less symetric
%
% http://bioinformatics.oxfordjournals.org/content/suppl/2014/11/10/btu747.DC1/Hi-Corrector_supplementary_materials.pdf
%
function [b,c] = Normalize(a)
	size_a = size(a,1);
	b = ones(1,size_a); %Row vector

	%Clean annoying stuff
	[f,xi] = ksdensity(sum(a),0:500:90000,'function','pdf');
	[pks,pksInd] = findpeaks(-f);
	if numel(pks) ~= 1
		fprintf('Warning, more than one cutoff! ');
		figure;plot(xi,f); % check this
		cutoff = xi(pksInd(1));
	else
		cutoff = xi(pksInd);
	end
	fprintf('Cutoff at %d... ',cutoff);
	t = sum(a);
	a(t<cutoff,:) = 0;
	a(:,t<cutoff) = 0;

	%Iterate. No stop condition? Odd.
	for i = 1:1000
		t = sum(a); %Row vector
		t(t==0) = Inf;
		a = a ./ (t'*t);
		b = b .* t;
	end
end
