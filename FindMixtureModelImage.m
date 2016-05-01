% This is for the paper. Meant to find nice looking diagonal 
function [a,a_tad,a_bg,a_cross] = FindMixtureModelImage(chrNum,prefix,res,win,diagNums)

dirPath = '/cs/grad/gilr02/cbio/htad-chain/output/';
fMatrix = [dirPath prefix sprintf('.matrix.chr%d.txt',chrNum)];
fDomains = [dirPath prefix sprintf('.domains.chr%d.txt',chrNum)];
	
blush = [0 66 123; 4 85 139; 10 97 154; 30 115 172; 40 130 200; 60 150 210]/255;
redsh = [171 55 54; 178 61 60; 196 66 66; 217 50 31; 231 79 78; 233 104 104]/255;
yelsh = [234 193 23; 242 187 102 ; 251 185 23 ; 251 177 23 ; 255 166 47 ; 233 171 23]/255;

if ischar(win)
	win = str2num(win);
end;
if ischar(res)
	res = str2num(res);
end;

MIN_DIAG = 1;%0 is the middle one; %SETTING MINDIAG TO HIGHER IS DANGEROUS
MAX_DIAG = win/res; %So that 2mb in 20k window will result in such and so

%Load
Log('Loading');
a = load(fMatrix);
d = load(fDomains);
Log();

a_size = size(a,2);

%Trim diagonals
Log('Creating diag matrix');
clear_diag = CreateDiagMatrix(a_size,MIN_DIAG,MAX_DIAG);
a_diag = a .* clear_diag;
Log();

%Use public TAD data to extract intraTAD stuff
Log('Extracting TAD/Bg matrix-masks'); %STORE THIS!!!! TAKES FOREVER
[a_bgmask,a_tadmask,a_crossmask] = GetTADsMask(d,MIN_DIAG,MAX_DIAG,a_size,res);
Log();

% Depracated
%a_bgmask = a_bgmask + a_crossmask; % A convetsion so that all non-basic TAD things will be bg
%a_crossmask = a_crossmask-a_crossmask; % Zero out a_cross so it will not disturb us
% Instead, we remove the crosses because they are a mixed state

a_tad = a_tadmask .* a_diag; %No tiles (xor), just TADs
a_tad(~a_tadmask)=NaN; 
a_bg = a_bgmask .* a_diag; %Background, no cross tads
a_bg(~a_bgmask)=NaN; 
a_cross = a_crossmask .* a_diag; %Cross TAD stuff
a_cross(~a_crossmask)=NaN; 

%To calculate priors, we set a_cross to 0 so that prior(tad) = #(tad)/#(tad)+#(bg)
[pr_tad, pr_bg, pr_crs] = GetPdTBC(a_tadmask,a_bgmask,a_crossmask.*0,MIN_DIAG,MAX_DIAG);%Model
Log();

figure('color','w'); %Used for multi diagonals
i = 1;
for diagNum = diagNums
	%All	
	[f,xi] = ksdensity( GetDiag(a,diagNum,1,1,1,1) , sort(GetDiag(a,diagNum,1,1,1,1)) , 'function','pdf' );
	%semilogx(exp(xi-1),f,'-','LineWidth',2.5,'Color',redsh(i,:)); hold on;

	%TAD
	diagTad = GetDiag(a_tad,diagNum,1,0,1,1);
	[f_t,xi_t] = ksdensity( diagTad , sort(diagTad) , 'function','pdf' );
	%semilogx(exp(xi_t-1),f_t*pr_tad(diagNum),'-','LineWidth',2.5,'Color',blush(i,:)); %This with priors
	semilogx(exp(xi_t-1),f_t,'-','LineWidth',2.5,'Color',blush(i,:)); %This with priors
	hold on;

	%Bg
	diagBg = GetDiag(a_bg,diagNum,1,1,1,1);
	[f_t,xi_t] = ksdensity( diagBg , sort(diagBg) , 'function','pdf' );
	%semilogx(exp(xi_t-1),f_t*pr_bg(diagNum),'-','LineWidth',2.5,'Color',yelsh(i,:));

	i = i + 1;
	xlabel('Interaction Intensity');
	ylabel('Probability Density');
	set(gca,'xtick',[1 10 100 1000],'xticklabel',[1 10 100 1000]);

	%figure;
	%qqplot(xi_t);
	% For saving lots of images: Was used for finding a good image
	%title(sprintf('%s chr%d diag %d res %d',prefix,chrNum,diagNum,res));
	%SaveFigure(gcf,sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/%s.%d.%d.%d.png',prefix,chrNum,diagNum,res));
	%SaveFigure(gcf,sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/%s.%d.%d.%d.fig',prefix,chrNum,diagNum,res));

	% Some statistical analisys
	%{
	%fprintf('Num\tType\tSW\tKS\tSKW\tKRT\r\n');
	data = GetDiag(a,diagNum,1,1,1,1);
	swo = swtest(data);
	kso = kstest(data);
	skw = skewness(data);
	krt = kurtosis(data);
	fprintf('%d\t%s\t%g\t%g\t%3.3g\t%3.3g\r\n',diagNum,'All',swo,kso,skw,krt);
	data = GetDiag(a_tad,diagNum,1,1,1,1);
	swo = swtest(data);
	kso = kstest(data);
	skw = skewness(data);
	krt = kurtosis(data);
	fprintf('%d\t%s\t%g\t%g\t%3.3g\t%3.3g\r\n',diagNum,'Tad',swo,kso,skw,krt);
	data = GetDiag(a_bg,diagNum,1,1,1,1);
	swo = swtest(data);
	kso = kstest(data);
	skw = skewness(data);
	krt = kurtosis(data);
	fprintf('%d\t%s\t%g\t%g\t%3.3g\t%3.3g\r\n',diagNum,'Bg',swo,kso,skw,krt);
	%}
end

return
%all this is not running lolz 

%Find gaussians
Log('Gaussianing');

a_gmm_tad = AnalyzeGMM(a_tad,MIN_DIAG,MAX_DIAG,1);
a_gmm_bg = AnalyzeGMM(a_bg,MIN_DIAG,MAX_DIAG,1);
meansTadDxn = zeros(size(a_gmm_tad));
meansBgDxn = zeros(size(a_gmm_bg));
sigmaTadDxn = zeros(size(a_gmm_tad));
sigmaBgDxn = zeros(size(a_gmm_bg));

a_gmm_all = AnalyzeGMM(a_diag,MIN_DIAG,MAX_DIAG,1);
meansAll = zeros(size(a_gmm_all));
for i = 1:numel(a_gmm_tad)
	meansTadDxn(i) = a_gmm_tad{i}.mu(1);
	meansBgDxn(i) = a_gmm_bg{i}.mu(1);
	sigmaTadDxn(i) = sqrt(a_gmm_tad{i}.Sigma(1));
	sigmaBgDxn(i) = sqrt(a_gmm_bg{i}.Sigma(1));
	
	meansAll(i) = a_gmm_all{i}.mu(1);
end
Log();

if 1 %Print Mean & Sigma for Dixon
	x_vals = log((MIN_DIAG:MAX_DIAG)*res);
	figure('color','w');
	plot(x_vals,meansAll,'LineWidth',2);
	hold on;
	plot(x_vals,meansTadDxn,'LineWidth',2);
	plot(x_vals,meansBgDxn,'LineWidth',2);
	legend('All','Tad','Bg');

	%loglog(x_vals,meansTadDxn+sigmaTadDxn,'r--');
	%loglog(x_vals,meansTadDxn-sigmaTadDxn,'r--');
	%loglog(x_vals,meansBgDxn+sigmaBgDxn,'b--');
	%loglog(x_vals,meansBgDxn-sigmaBgDxn,'b--');
	title('Mean vs. Distance');
	xlabel('logDistance (bp)');
	ylabel('Mean');
	fprintf('Done\r\n');
end

%Plot PdT, pdB, pdC
if 0
	x_vals = (MIN_DIAG:MAX_DIAG)*res;
	figure; hold on;
	plot(x_vals,pr_tad,'LineWidth',2);
	plot(x_vals,pr_bg,'LineWidth',2);
	legend('TAD','Bg');
	xlabel('Distance (bp)');
	ylabel('Prob');
	title('Pr(TAD), Pr(Bg) vs. Distance');
end

end
