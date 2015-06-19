function [a_ret] = TAD(chrNumber)
DATA_FILE = sprintf('HiC-CSV-Matrices/IMR90_chr%d_40k', chrNumber);
MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = 75;%55
BLOCK_SIZE = 40000;

%Quack, initialize return value
a_ret = 0;

%Load
Log('Loading');
a=load(DATA_FILE);
a_size = size(a,1);
Log();

%Normalize
%% The problematic cells are killed by b, having inf in their location!
%% Any following code MUST drop zero value cells.
%% This also kills other zeros cells
Log('Normalizing');
norm_vec_filename = sprintf('HiC-CSV-Matrices/normalization_vector_IMR90_chr%d_40k.mat', chrNumber);
if exist(norm_vec_filename,'file')
	fprintf('Norm vector found');

	load(norm_vec_filename,'b');
else
	fprintf('Recalculating norm vector');
	b = Normalize(a);
	fprintf('Saving to dist');
	save(norm_vec_filename,'b');
end
a = a ./ (b'*b);
% ( sum(a)==0) ) == (b==Inf)
% These are the indices of interest
% ind = (b!=Inf)
% a_clean = a(ind,ind)
% a_orig(ind,ind) = a_clean
% When removing, TAD structure is curios as well...
Log();

if 0 %Normalize with evil unethical method
	Log('Equalizing Spectrum...');
	a_eq = EqualizeSpectrum(a,MIN_DIAG,MAX_DIAG);
	Log();	
end
%a_ret = a; return;

%Trim diagonal
Log('Creating diag matrix');
clear_diag = triu(tril(ones(size(a)),MAX_DIAG),MIN_DIAG);
a_diag = a .* clear_diag;
Log();

%Use public TAD data to extract intraTAD stuff
Log('Extracting TADs');
[a_bgmap,a_tadmap,a_crossmap] = GetTADs(chrNumber,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE);
a_tad = a_tadmap .* a_diag;%No tiles, just TADs
a_bound = a_bgmap .* a_diag; %For bounderies
a_cross = a_crossmap .* a_diag; %Cross TAD stuff
[pr_tad, pr_bg pr_crs] = GetPdTBC(a_diag,a_tad,a_bound,a_cross,MIN_DIAG,MAX_DIAG);%Model
Log();

%Find gaussians
Log('Gaussianing');
a_gmm_tad = AnalyzeGMM(a_tad,MIN_DIAG,MAX_DIAG,1);
a_gmm_bound = AnalyzeGMM(a_bound,MIN_DIAG,MAX_DIAG,1);
a_gmm_cross = AnalyzeGMM(a_cross,MIN_DIAG,MAX_DIAG,1);
meansTadDxn = zeros(size(a_gmm_tad));
meansBgDxn = zeros(size(a_gmm_bound));
meansCrsDxn = zeros(size(a_gmm_cross));
sigmaTadDxn = zeros(size(a_gmm_tad));
sigmaBgDxn = zeros(size(a_gmm_bound));
sigmaCrsDxn = zeros(size(a_gmm_cross));
for i = 1:numel(a_gmm_tad)
	meansTadDxn(i) = a_gmm_tad{i}.mu(1);
	meansBgDxn(i) = a_gmm_bound{i}.mu(1);
	meansCrsDxn(i) = a_gmm_cross{i}.mu(1);
	sigmaTadDxn(i) = a_gmm_tad{i}.Sigma(1);
	sigmaBgDxn(i) = a_gmm_bound{i}.Sigma(1);
	sigmaCrsDxn(i) = a_gmm_cross{i}.Sigma(1);
end
Log();

if 0 %Print Mean & Sigma for Dixon
x_vals = (MIN_DIAG:MAX_DIAG)*BLOCK_SIZE;
figure;
hold on;
loglog(x_vals,meansTadDxn,'rx-','LineWidth',2);
loglog(x_vals,meansBgDxn,'bx-','LineWidth',2);
loglog(x_vals,meansCrsDxn,'gx-','LineWidth',2);
legend('Tad','Bound','Cross');
loglog(x_vals,meansTadDxn+sigmaTadDxn,'r--');
loglog(x_vals,meansTadDxn-sigmaTadDxn,'r--');
loglog(x_vals,meansBgDxn+sigmaBgDxn,'b--');
loglog(x_vals,meansBgDxn-sigmaBgDxn,'b--');
loglog(x_vals,meansCrsDxn+sigmaCrsDxn,'g--');
loglog(x_vals,meansCrsDxn-sigmaCrsDxn,'g--');
title('Mean vs. Distance');
xlabel('Distance (bp)');
ylabel('Mean');
fprintf('Done\r\n');
end

%Generate GMM data
if 0
	TryTripleGMM(a_diag,MIN_DIAG,MAX_DIAG);
end

%Generate & Parse 2 GMM data
if 0
Log('GMMing (2 Param)');
a_gmm = AnalyzeGMM(a_diag,MIN_DIAG,MAX_DIAG,2); %CACH THIS WOWZER LIKE NORMALIZE DWAG YEH
%I mean cache the results of this long calculation
[meansBg, meansIn, sigmaBg, sigmaIn, probBg, probIn] = ParseGM(a_gmm, true);
Log();
end

%Get P_D_T(N) and P_D_B(N)
if 1
%Smooth data before continuing
%meansBg = smooth(meansBg);
%meansIn = smooth(meansIn);
%sigmaBg = smooth(sigmaBg);
%sigmaIn = smooth(sigmaIn);

%box = 1300:1800; %chr5
box = 1100:1700; %chr3
box = 400:a_size;

[a_llr,a_pdt,a_pdb] = PosteriorHeatmap(a_diag,MIN_DIAG,MAX_DIAG,pr_bg,meansBgDxn,sigmaBgDxn,pr_tad,meansTadDxn,sigmaTadDxn);
%[a_llr,a_pdt,a_pdb] = PosteriorHeatmap(a_diag,MIN_DIAG,MAX_DIAG,probBg,meansBg,sigmaBg,probIn,meansIn,sigmaIn);

%DisplayHeatmap( triu(a_pdt) + tril(a_pdb') , [0 1] , box );
%title(sprintf('Pd(TAD) & Pd(BG), chr%d',chrNumber));

DisplayHeatmap( triu(a_llr) + tril(a_tadmap' + a_crossmap') , [-3 3], box);
title(sprintf('Logaritmic liklyhood ratio, chr%d',chrNumber));

Log('Logaritmic Liklyhooding');
[a_pyrsky, a_pyr, a_sky] = PyramidSky(a_llr,-a_llr,MAX_DIAG,2);
Log();

Log('Deriving TADs from matrix');
[a_t_mp, a_t_sl] = DynProgTAD(a_pyrsky,box(1),box(end),a_pyr);
a_ret = zeros(size(box,2),size(box,2),2);
a_ret(:,:,1) = a_t_mp;
a_ret(:,:,2) = a_t_sl;
Log();

a_log_diag = log(a_diag+1);
a_log_diag = max(a_pyrsky(:))*(a_log_diag-min(a_log_diag(:)))/(max(a_log_diag(:))-min(a_log_diag(:)))+min(a_pyrsky(:));
%DisplayHeatmap(triu(a_pyrsky) + tril(a_log_diag'),[min(a_pyrsky(:)) max(a_pyrsky(:))], box);
title(sprintf('LLR vs. log(heatmap), chr%d',chrNumber));

%Here can do also -
% rng = [453, 509];
% asum = rng(1):rng(2); for i=rng(1):rng(2) vl=a(rng(1),i)+a(i+1,rng(2)); asum(i-rng(1)+1)=vl; end; end;
% figure; plot(asum); title(sprintf('%d-%d, chr%d',rng(1),rng(2),chrNumber));

%Also, add cost parameter to DynProgTAD. Then find tad solutions for changing cost values

end

%Plot PdT, pdB, pdC
if 0
x_vals = (MIN_DIAG:MAX_DIAG)*BLOCK_SIZE;
figure; hold on;
plot(x_vals,pr_tad);
plot(x_vals,pr_bg);
plot(x_vals,pr_crs);
legend('TAD','Bg','Cross');
xlabel('Distance (bp)');
ylabel('Prob');
title(sprintf('Pr(TAD), Pr(Bg), Pr(Cross) vs. Distance, for chr %d',chrNumber));
end

%Extract probability for TAD diagonals only
if 0
figure;
iplt = 1;
for i = [2,5,10,20,30,45]
	subplot(2,3,iplt);
	hold on;

	dgn_clean = GetDiag(a_tad,i,1,1,1,0);
	if numel(dgn_clean) > 0
		[f,xi] = ksdensity(dgn_clean,dgn_clean,'function','pdf');
		h = qqplot(dgn_clean);
		set(h(3),'linewidth',1,'color',[0 0 1]);
		set(h(2),'linewidth',1,'color',[0 0 1]);
		set(h(1),'marker','^','markersize',8,'markeredgecolor',[0 0 1]);
		%scatter(xi,f,11,'filled');
	end

	dgn_clean = GetDiag(a_bound,i,1,1,1,0);
	if numel(dgn_clean) > 0
		[f,xi] = ksdensity(dgn_clean,dgn_clean,'function','pdf');
		h = qqplot(dgn_clean);
		set(h(3),'linewidth',1,'color',[0 1 0]);
		set(h(2),'linewidth',1,'color',[0 1 0]);
		set(h(1),'marker','^','markersize',8,'markeredgecolor',[0 1 0]);
		%scatter(xi,f,11,'filled');
	end

	dgn_clean = GetDiag(a_cross,i,1,1,1,0);
	if numel(dgn_clean) > 0
		[f,xi] = ksdensity(dgn_clean,dgn_clean,'function','pdf');
		h = qqplot(dgn_clean);
		set(h(3),'linewidth',1,'color',[1 0 0]);
		set(h(2),'linewidth',1,'color',[1 0 0]);
		set(h(1),'marker','^','markersize',8,'markeredgecolor',[1 0 0]);
		%scatter(xi,f,11,'filled');
	end

	title(sprintf('Distance=%d (i=%d)',i*BLOCK_SIZE,i));
	iplt = iplt+1;
end
hb = plot(nan,nan,'.b');
hg = plot(nan,nan,'.g');
hr = plot(nan,nan,'.r');
legend([hb hg hr],'tad','bg','cross');
suptitle('Probability density of interactions strength, for several diagonals');
end

end
