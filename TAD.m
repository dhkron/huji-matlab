function [a_ret] = TAD(chrNumber,res)
MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = floor(65*40000/res);%55
BLOCK_SIZE = res;%40k PLEASE

%Quack, initialize return value
a_ret = 0;

%Load & Normalize
Log('Loading & Normalizng');
if res==40000
	[~,a] = RegularLoader(chrNumber);
	%[~,a] = LizLoader(chrNumber);
	a = a/2;
else
	[~,a] = LibermanLoader(chrNumber,res);
	a = a;
end
a_size = size(a,1);
Log();

if 0 %Normalize with evil unethical method
	Log('Equalizing Spectrum...');
	a_eq = EqualizeSpectrum(a,MIN_DIAG,MAX_DIAG);
	Log();	
end
%a_ret = a; return;

%Trim diagonal
Log('Creating diag matrix');
clear_diag = CreateDiagMatrix(a_size,MIN_DIAG,MAX_DIAG);
a_diag = a .* clear_diag;
Log();

if 0
	%Use public TAD data to extract intraTAD stuff
	Log('Extracting TADs'); %STORE THIS!!!! TAKES FOREVER
	[a_bgmap,a_tadmap,a_crossmap] = GetTADs(chrNumber,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE);
	a_tad = a_tadmap .* a_diag;%No tiles, just TADs
	a_bound = a_bgmap .* a_diag; %For bounderies
	a_cross = a_crossmap .* a_diag; %Cross TAD stuff
	[pr_tad, pr_bg pr_crs] = GetPdTBC(a_tadmap,a_bgmap,a_cross.*0,MIN_DIAG,MAX_DIAG);%Model
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
		sigmaTadDxn(i) = sqrt(a_gmm_tad{i}.Sigma(1));
		sigmaBgDxn(i) = sqrt(a_gmm_bound{i}.Sigma(1));
		sigmaCrsDxn(i) = sqrt(a_gmm_cross{i}.Sigma(1));
	end
	Log();
else
	t = load('Median40k.65.mat');
	meansTadDxn = t.S.mt;
	meansBgDxn = t.S.mb;
	sigmaTadDxn = t.S.st;
	sigmaBgDxn = t.S.sb;
	%Here interpolate
	%interp1(x,y,nx)
	%x = in basepairs and stuff
	%y = meansTadDxn and stuff
	%nx = basepairs
	%problem - need to extrapolate at diag=1
	x_vals = (1:65)*40000;
	x_vals_new = (MIN_DIAG:MAX_DIAG)*res;
	meansTadDxn = interp1(x_vals,meansTadDxn,x_vals_new);
	meansBgDxn = interp1(x_vals,meansBgDxn,x_vals_new);
	sigmaTadDxn = interp1(x_vals,sigmaTadDxn,x_vals_new);
	sigmaBgDxn = interp1(x_vals,sigmaBgDxn,x_vals_new);
end

if 0 %Print Mean & Sigma for Dixon
x_vals = (MIN_DIAG:MAX_DIAG)*BLOCK_SIZE;
figure;
loglog(x_vals,meansTadDxn,'rx-','LineWidth',2);
hold on;
loglog(x_vals,meansBgDxn,'bx-','LineWidth',2);
%loglog(x_vals,meansCrsDxn,'gx-','LineWidth',2);
legend('Tad','Bound');%,'Cross');
loglog(x_vals,meansTadDxn+sigmaTadDxn,'r--');
loglog(x_vals,meansTadDxn-sigmaTadDxn,'r--');
loglog(x_vals,meansBgDxn+sigmaBgDxn,'b--');
loglog(x_vals,meansBgDxn-sigmaBgDxn,'b--');
%loglog(x_vals,meansCrsDxn+sigmaCrsDxn,'g--');
%loglog(x_vals,meansCrsDxn-sigmaCrsDxn,'g--');
title('Mean vs. Distance');
xlabel('Distance (bp)');
ylabel('Mean');
fprintf('Done\r\n');
end

%Plot PdT, pdB, pdC
if 0
x_vals = (MIN_DIAG:MAX_DIAG)*BLOCK_SIZE;
figure; hold on;
plot(x_vals,pr_tad,'LineWidth',2);
plot(x_vals,pr_bg,'LineWidth',2);
%plot(x_vals,pr_crs);
legend('TAD','Bg');
%legend('TAD','Bg','Cross');
xlabel('Distance (bp)');
ylabel('Prob');
title(sprintf('Pr(TAD), Pr(Bg) vs. Distance, for chr %d',chrNumber));
%title(sprintf('Pr(TAD), Pr(Bg), Pr(Cross) vs. Distance, for chr %d',chrNumber));
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
[meansBg, meansIn, sigmaBg, sigmaIn, pr_bg, pr_td] = ParseGM(a_gmm, true);
Log();
end


%Smooth data before continuing
if 1
	%meansBg = smooth(meansBg);
	%meansIn = smooth(meansIn);
	%sigmaBg = smooth(sigmaBg);
	%sigmaIn = smooth(sigmaIn);
end

%box = 1300:1800; %chr5
%box = 1100:1700; %chr3
%box = floor(400*40000/BLOCK_SIZE):a_size;
%box = 1:a_size;
%box = 500:800;%chr19 liz
box = 200:500;%chr19 liz

%Posterior or likelihood
if 0
	[a_lpr,a_pdt,a_pdb] = PosteriorHeatmap(a_diag,MIN_DIAG,MAX_DIAG,pr_bg,meansBgDxn,sigmaBgDxn,pr_tad,meansTadDxn,sigmaTadDxn);
else
	[a_lpr,a_pdt,a_pdb] = LiklihoodHeatmap(a_diag,MIN_DIAG,MAX_DIAG,meansBgDxn,sigmaBgDxn,meansTadDxn,sigmaTadDxn);
end

%DisplayHeatmap( triu(a_pdt) + tril(a_pdb') , [0 1] , box );
%title(sprintf('Pd(TAD) & Pd(BG), chr%s',chrNumber));

%DisplayHeatmap( triu(log(a_diag+1)) + tril(a_tadmap' + a_crossmap') , [-4 4], box);

Log('Logaritmic Liklyhooding');
a_lpr(a_lpr== inf) = 0;
a_lpr(a_lpr==-inf) = 0;
[a_pyrsky, a_pyr, a_sky] = PyramidSky(a_lpr,-a_lpr,MAX_DIAG,2);
Log();

%figure; hold on;
%rng = [740 816];%Single
%rng = [453 509];%Good 2
%rng = [235 301]; %should be 2
%a_pyrsky_box = a_pyrsky(box,box);
%tadvalue = a_pyrsky_box(rng(1) , rng(2) );
%plot([1 rng(2)-rng(1)+1] , [tadvalue tadvalue],'--');

a_old_mask = ones(1,a_size);
for cost = [0 100 300 1000 3000 10000]
	%DisplayHeatmap( triu(a_lpr) + tril(a_lpr') , [-4 4], box);
	%f = DisplayHeatmap(log(a), 0, box, 'red');
	f = DisplayHeatmap(a_lpr, [-2 4], box); %lpr = llr
	title(sprintf('LogLikelihood + TAD mapping, chr%s cost=%d res=%d',chrNumber,cost,BLOCK_SIZE));
	Log(sprintf('Deriving TADs from matrix, cost = %d',cost));
	%[a_t_mp, a_t_sl, a_bnd_mask] = DynProgTAD(a_pyrsky,box(1),box(end),cost,a_old_mask);
	[a_t_mp, a_t_sl, a_bnd_mask] = DynProgTAD(a_pyrsky,box(1),box(end),cost,a_old_mask);
	a_old_mask = a_bnd_mask;
	a_ret = zeros(size(box,2),size(box,2),2);
	a_ret(:,:,1) = a_t_mp;
	a_ret(:,:,2) = a_t_sl;
	%caxis([-9 -5]);
	caxis([-1 1]);
	if 1
		figure_save_path=sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/Chr%d_%d_%d_cost_%d_llr.png',chrNumber,box(1),box(end),cost);
		SaveFigure(f,figure_save_path);
	end
	Log();

	%Here can do also -
%	asum = rng(1):rng(2)-1;
%	for i=rng(1):rng(2)-1
%		vl=a_t_mp(rng(1),i)+a_t_mp(i+1,rng(2))-cost;
%		asum(i-rng(1)+1)=vl;
%	end;
%	plot(asum); title(sprintf('%d-%d, chr%s',rng(1),rng(2),chrNumber));
end
%legend('quack','100000' ,'10000' ,'3000' ,'1000' ,'300' ,'10' ,'0');

%Display supersum vs. heatmap
%a_log_diag = log(a_diag+1);
%a_log_diag = max(a_pyrsky(:))*(a_log_diag-min(a_log_diag(:)))/(max(a_log_diag(:))-min(a_log_diag(:)))+min(a_pyrsky(:));
%DisplayHeatmap(triu(a_pyrsky) + tril(a_log_diag'),[min(a_pyrsky(:)) max(a_pyrsky(:))], box);
%title(sprintf('Supersum vs. log(heatmap), chr%s',chrNumber));

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
