% This is for the paper. Meant to find nice looking diagonal 
function [a,a_tad,a_bg,a_cross] = FindMixtureModelImage(chrNum,prefix,res,win,diagNums)

dirPath = '/cs/grad/gilr02/cbio/htad-chain/output/';
fMatrix = [dirPath prefix sprintf('.matrix.chr%d.txt',chrNum)];
fDomains = [dirPath prefix sprintf('.domains.chr%d.txt',chrNum)];

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
a_bg = a_bgmask .* a_diag; %Background, no cross tads
a_cross = a_crossmask .* a_diag; %Cross TAD stuff

%To calculate priors, we set a_cross to 0 so that prior(tad) = #(tad)/#(tad)+#(bg)
[pr_tad, pr_bg, pr_crs] = GetPdTBC(a_tadmask,a_bgmask,a_crossmask.*0,MIN_DIAG,MAX_DIAG);%Model
Log();

%fprintf('Num\tType\tSW\tKS\tSKW\tKRT\r\n');
for diagNum = diagNums
	%{
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

	figure;
	title(sprintf('%s chr%d diag %d res %d',prefix,chrNum,diagNum,res));
	hold on;
	[f,xi] = ksdensity( GetDiag(a,diagNum,1,1,1,1) , sort(GetDiag(a,diagNum,1,1,1,1)) , 'function','pdf' );
	plot(xi,f,'-');
	[f_t,xi_t] = ksdensity( GetDiag(a_tad,diagNum,1,1,1,1) , sort(GetDiag(a_tad,diagNum,1,1,1,1)) , 'function','pdf' );
	plot(xi_t,f_t*pr_tad(diagNum),'-');
	[f_b,xi_b] = ksdensity( GetDiag(a_bg,diagNum,1,1,1,1) , sort(GetDiag(a_bg,diagNum,1,1,1,1)) , 'function','pdf' );
	plot(xi_b,f_b*pr_bg(diagNum),'-');
	SaveFigure(gcf,sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/%s.%d.%d.%d.png',prefix,chrNum,diagNum,res));
	SaveFigure(gcf,sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/%s.%d.%d.%d.fig',prefix,chrNum,diagNum,res));
end

return
%all this is not running lolz 
%{
%Find gaussians
Log('Gaussianing');

a_gmm_tad = AnalyzeGMM(a_tad,MIN_DIAG,MAX_DIAG,1);
a_gmm_bg = AnalyzeGMM(a_bg,MIN_DIAG,MAX_DIAG,1);
meansTadDxn = zeros(size(a_gmm_tad));
meansBgDxn = zeros(size(a_gmm_bg));
sigmaTadDxn = zeros(size(a_gmm_tad));
sigmaBgDxn = zeros(size(a_gmm_bg));
for i = 1:numel(a_gmm_tad)
	meansTadDxn(i) = a_gmm_tad{i}.mu(1);
	meansBgDxn(i) = a_gmm_bg{i}.mu(1);
	sigmaTadDxn(i) = sqrt(a_gmm_tad{i}.Sigma(1));
	sigmaBgDxn(i) = sqrt(a_gmm_bg{i}.Sigma(1));
end
Log();

Log('Saving');
save(fOut,'meansTadDxn','meansBgDxn','sigmaTadDxn','sigmaBgDxn','pr_tad','pr_bg');
Log();

if 0 %Print Mean & Sigma for Dixon
	x_vals = (MIN_DIAG:MAX_DIAG)*res;
	figure;
	loglog(x_vals,meansTadDxn,'rx-','LineWidth',2);
	hold on;
	loglog(x_vals,meansBgDxn,'bx-','LineWidth',2);
	legend('Tad','Bound');

	loglog(x_vals,meansTadDxn+sigmaTadDxn,'r--');
	loglog(x_vals,meansTadDxn-sigmaTadDxn,'r--');
	loglog(x_vals,meansBgDxn+sigmaBgDxn,'b--');
	loglog(x_vals,meansBgDxn-sigmaBgDxn,'b--');
	title('Mean vs. Distance');
	xlabel('Distance (bp)');
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
%}
