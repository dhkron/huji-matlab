function [meansTadDxn,sigmaTadDxn,meansBgDxn,sigmaBgDxn,meansDiagDxn,sigmaDiagDxn] = UniversalMeanRule(chrs)
res=40000;
MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = floor(65*40000/res);%55
BLOCK_SIZE = res;%40k PLEASE

meansTadDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);
meansBgDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);
meansDiagDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);
sigmaTadDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);
sigmaBgDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);
sigmaDiagDxn = zeros(MAX_DIAG-MIN_DIAG+1,22);

figure;
%Loopy
for chrNumber=chrs
	Log(sprintf('Chr %d',chrNumber));
	Log();

	%Quack, initialize return value
	a_ret = 0;

	%Load & Normalize
	Log('Loading & Normalizng');
	[~,a] = RegularLoader(chrNumber);
	a_size = size(a,1);
	Log();

	%Trim diagonal
	Log('Creating diag matrix');
	clear_diag = CreateDiagMatrix(a_size,MIN_DIAG,MAX_DIAG);
	a_diag = a .* clear_diag;
	Log();

	%Use public TAD data to extract intraTAD stuff
	Log('Extracting TADs'); %STORE THIS!!!! TAKES FOREVER
	[a_bgmap,a_tadmap,~] = GetTADs(chrNumber,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE);
	a_tad = a_tadmap .* a_diag;%No tiles, just TADs
	a_bound = a_bgmap .* a_diag; %For bounderies
	Log();

	%Find gaussians
	Log('Gaussianing');
	a_gmm_tad = AnalyzeGMM(a_tad,MIN_DIAG,MAX_DIAG,1);
	a_gmm_bound = AnalyzeGMM(a_bound,MIN_DIAG,MAX_DIAG,1);
	a_gmm_diag = AnalyzeGMM(a_diag,MIN_DIAG,MAX_DIAG,1);
	for i = 1:numel(a_gmm_tad)
		meansTadDxn(i,chrNumber) = a_gmm_tad{i}.mu(1);
		meansBgDxn(i,chrNumber) = a_gmm_bound{i}.mu(1);
		meansDiagDxn(i,chrNumber) = a_gmm_diag{i}.mu(1);
		sigmaTadDxn(i,chrNumber) = sqrt(a_gmm_tad{i}.Sigma(1));
		sigmaBgDxn(i,chrNumber) = sqrt(a_gmm_bound{i}.Sigma(1));
		sigmaDiagDxn(i,chrNumber) = sqrt(a_gmm_diag{i}.Sigma(1));
	end
	Log();

	if 0 %Display tad,bg,diag with sigma margin
		x_vals = (MIN_DIAG:MAX_DIAG)*BLOCK_SIZE;
		loglog(x_vals,meansTadDxn(:,chrNumber),'rx-','LineWidth',2);
		hold on;
		loglog(x_vals,meansBgDxn(:,chrNumber),'bx-','LineWidth',2);
		loglog(x_vals,meansDiagDxn(:,chrNumber),'gx-','LineWidth',2);
		legend('Tad','Bound');
		loglog(x_vals,meansTadDxn(:,chrNumber)+sigmaTadDxn(:,chrNumber),'r--');
		loglog(x_vals,meansTadDxn(:,chrNumber)-sigmaTadDxn(:,chrNumber),'r--');
		loglog(x_vals,meansBgDxn(:,chrNumber)+sigmaBgDxn(:,chrNumber),'b--');
		loglog(x_vals,meansBgDxn(:,chrNumber)-sigmaBgDxn(:,chrNumber),'b--');
		loglog(x_vals,meansDiagDxn(:,chrNumber)+sigmaDiagDxn(:,chrNumber),'g--');
		loglog(x_vals,meansDiagDxn(:,chrNumber)-sigmaDiagDxn(:,chrNumber),'g--');
		title('Mean vs. Distance');
		xlabel('Distance (bp)');
		ylabel('Mean');
		fprintf('Done\r\n');
	end

	if 1 %Plot mean vs. diag
		loglog(meansDiagDxn(:,chrNumber),meansTadDxn(:,chrNumber),'rx-','LineWidth',2);
		hold on;
		loglog(meansDiagDxn(:,chrNumber),meansBgDxn(:,chrNumber),'bx-','LineWidth',2);
		legend('Tad','Background');
		title('Tad/Background Mean vs. Diag Mean intensity');
		xlabel('Diagonal Mean Intensity (log scale)');
		ylabel('TAD/Background Mean Intensity (log scale)');
		fprintf('Done\r\n');
	end
end

return

figure;
meansTadMed = median(meansTadDxn,2);
meansBgMed = median(meansBgDxn,2);
sigmaTadMed = median(sigmaTadDxn,2);
sigmaBgMed = median(sigmaBgDxn,2);
loglog(x_vals,meansTadMed,'rx-','LineWidth',2);
hold on;
loglog(x_vals,meansBgMed,'bx-','LineWidth',2);
loglog(x_vals,meansTadMed+sigmaTadMed,'r--');
loglog(x_vals,meansTadMed-sigmaTadMed,'r--');
loglog(x_vals,meansBgMed+sigmaBgMed,'b--');
loglog(x_vals,meansBgMed-sigmaBgMed,'b--');

return

end
