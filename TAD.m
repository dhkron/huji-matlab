function [a_ret] = TAD(chrNumber)
DATA_FILE = sprintf('HiC-CSV-Matrices/IMR90_chr%d_40k', chrNumber);
MIN_DIAG = 1;%0 is the middle one
MAX_DIAG = 55;%55
BLOCK_SIZE = 40000;

%Quack, initialize return value
a_ret = 0;

%For directionality index
DI_TRIM = 1/250.0;%Both work, interestingly. 1e-10%1/250.0;
TRANS_GUESS = [
[0.9000 0.0900 0.0100];
[0.1000 0.8000 0.1000];
[0.0100 0.0900 0.9000];
];
EMIS_GUESS = [
[0.7130 0.2754 0.0058 0.0058 0];
[0 0 0.9828 0.0172 0];
[0 0.0055 0.0247 0.3205 0.6493];
];

%Load
fprintf('Loading... ');
a=load(DATA_FILE);
fprintf('Done\r\n');

%Normalize
%% The problematic cells are killed by b, having inf in their location!
%% Any following code MUST drop zero value cells.
%% This also kills other zeros cells
fprintf('Normalizing... ');
norm_vec_filename = sprintf('HiC-CSV-Matrices/normalization_vector_IMR90_chr%d_40k.mat', chrNumber);
if exist(norm_vec_filename,'file')
	fprintf('Norm vector found... ');
	load(norm_vec_filename,'b');
else
	fprintf('Recalculating norm vector... ');
	b = Normalize(a);
	fprintf('Saving to dist... ');
	save(norm_vec_filename,'b');
end
a = a ./ (b'*b);
fprintf('Done\r\n');

%Sort diagonals to see what happens
figure;
for i = MIN_DIAG+5%MAX_DIAG 
	dgn = diag(a,i);
	dgn = dgn(dgn~=0);
	dgn = sort(dgn);
	dgn = log(dgn);
	plot(dgn);
	hold on;
end
hold off;
a_ret = a; return;

%Trim diagonal
fprintf('Creating diag matrix... ');
clear_diag = triu(tril(ones(size(a)),MAX_DIAG),MIN_DIAG);
a_diag = a .* clear_diag;
fprintf('Done\r\n');

%Generate GMM data
fprintf('GMMing... ');
a_size = size(a,1);
a_gmm = {};
a_pdt = zeros(a_size);
a_pdb = zeros(a_size);
warning('off','stats:gmdistribution:FailedToConverge');
for i = MIN_DIAG:MAX_DIAG
	dgn = diag(a_diag,i);
	dgn_clean = dgn(dgn~=0); %It is important to fit WITHOUT the zeros
	dgn_clean = log(dgn_clean+1);

	options = statset('MaxIter',1000);
	GMModel = fitgmdist(dgn_clean,2,'Options',options);
	a_gmm{i}=GMModel;

	if GMModel.mu(1) < GMModel.mu(2)
		i_bg = 1;
		i_in = 2;
	else
		i_bg = 2;
		i_in = 1;
	end
	start_pos = a_size*i + 1;
	diag_elem = log(diag(a_diag,i)+1); %Still maps zero to zero

	pdt_diag = normpdf(diag_elem,GMModel.mu(i_in),GMModel.Sigma(i_in))*GMModel.ComponentProportion(i_in);
	pdt_diag(diag_elem==0) = 0;
	a_pdt(start_pos:a_size+1:end) = pdt_diag;

	pdb_diag = normpdf(diag_elem,GMModel.mu(i_bg),GMModel.Sigma(i_bg))*GMModel.ComponentProportion(i_bg);
	pdb_diag(diag_elem==0) = 0;
	a_pdb(start_pos:a_size+1:end) = pdb_diag;

	[~,msgid] = lastwarn;
	if strcmp(msgid, 'stats:gmdistribution:FailedToConverge')
		fprintf('Failed to converge at diag %d\r\n',i);
		warning('')
	end
end
w = warning('query','last');
fprintf('Done\r\n');

a_unified = log(a_pdt)-log(a_pdb);

%Use public TAD data to extract intraTAD stuff
fprintf('Extracting TADs... ');
domains = load(sprintf('TADs/domains.IMR90.chr%d', chrNumber));
domains = floor(domains/BLOCK_SIZE)+1;
domains_map = zeros(a_size);
domains_map_xor = zeros(a_size);
domains_map_one = zeros(a_size);
%For calculating TAD strength vs. Length
tad_size = zeros(1,a_size);
tad_color = zeros(1,a_size);
tad_index = 1;
for i = MIN_DIAG:MAX_DIAG
	start_pos = a_size*i + 1;
	domains_diag_xor = zeros(1,a_size-i);
	domains_diag_one = zeros(1,a_size-i);
	for domain = domains'
		domain(2) = min(domain(2),a_size); %Fix overflow

		bnd_start = domain(1);
		bnd_end = domain(2)-i;
		min_bnd = min(bnd_start,bnd_end);
		min_bnd = max(1, min_bnd);
		max_bnd = max(bnd_start,bnd_end);
		max_bnd = min(max_bnd, a_size-i);
		domains_diag_xor(min_bnd:max_bnd) = ~domains_diag_xor(min_bnd:max_bnd);

		min_bnd = max(1, bnd_start);
		max_bnd = min(bnd_end, a_size-i);
		domains_diag_one(min_bnd:max_bnd) = 1; %Do not extend TADs to tiles

		%Calcualte mean color
		tad_sum = sum(sum((a_unified(domain(1):domain(2),domain(1):domain(2))))); %Can use log here
		tad_count = numel(a_unified(domain(1):domain(2),domain(1):domain(2)));
		tad_size(tad_index) = domain(2)-domain(1);
		tad_color(tad_index) = tad_sum/tad_count;
		tad_index = tad_index + 1;
	end
	domains_map_xor(start_pos:a_size+1:end) = domains_diag_xor;
	domains_map_one(start_pos:a_size+1:end) = domains_diag_one;
	%domains_map(start_pos:a_size+1:end) = round(rand(1,a_size-i)); %Random model
end
domains_map_crs = domains_map_xor-domains_map_one; %Just TAD interactions, not TADs

%Plot TAD strength vs. Length
if 1 
figure;
scatter(tad_size*BLOCK_SIZE,tad_color);
title('Mean interaction strength vs. TAD length');
xlabel('Length (bp)');
ylabel('Mean Interaction strength');
set(gca,'xscale','log');
set(gca,'yscale','log');
coeffs = polyfit(log(tad_size*BLOCK_SIZE),log(tad_color),1);
hold on;
plot(tad_size*BLOCK_SIZE,((tad_size*BLOCK_SIZE).^coeffs(1))*exp(coeffs(2)));
hold off;
legend('Original Data',sprintf('Line Fit ax^b\r\na=%d\r\nb=%d',exp(coeffs(2)),coeffs(1)));
end

%Select the proper separation model
domains_map = domains_map_one; %No tiles, just TADs
bounds_map = 1-domains_map_xor; %Whatever is not tads or interactions
a_tad = domains_map .* a_diag;
a_bound = bounds_map .* a_diag; %For bounderies
fprintf('Done\r\n');

%Extract probability for TAD diagonals only
if 0
figure;
iplt = 1;
for i = [2,5,10,20,30,45]
	dgn = diag(a_tad,i);
	dgn_clean = dgn(dgn~=0);
	dgn_clean = log(dgn_clean+1);
	[f,xi] = ksdensity(dgn_clean,dgn_clean,'function','pdf');

	subplot(2,3,iplt);
	scatter(xi,f,11,'filled');
	hold on;

	dgn = diag(a_bound,i);
	dgn_clean = dgn(dgn~=0);
	dgn_clean = log(dgn_clean+1);
	[f,xi] = ksdensity(dgn_clean,dgn_clean,'function','pdf');

	scatter(xi,f,11,'filled');
	title(sprintf('Distance=%d (i=%d)',i*BLOCK_SIZE,i));
	legend('intra-domains','bounderies');
	iplt = iplt+1;
end
suptitle('Probability density of interactions strength, for several diagonals');
end

%QQPlots for several diags
if 0
dgn = diag(a_bound,2);
dgn = dgn(dgn~=0);
figure;qqplot(log(dgn+1));title('qq plot, bound, 80k');
dgn = diag(a_tad,2);
dgn = dgn(dgn~=0);
figure;qqplot(log(dgn+1));title('qq plot, tad, 80k');
dgn = diag(a_tad,20);
dgn = dgn(dgn~=0);
figure;qqplot(log(dgn+1));title('qq plot, tad, 800k');
end

%Plot Public TAD mapping & Original image
dispbox = 1:a_size;
if 1
a_diag_disp = triu(log(a_diag))+tril(domains_map_one' + domains_map_xor');
a_diag_disp = a_diag_disp(dispbox,dispbox);
figure; imagesc(a_diag_disp); axis equal; colorbar; caxis([0 8]);
title('log(Heatmap) vs. Dixon TAD Mapping');
end

if 1
a_unif_disp = triu(a_unified)+tril(domains_map_one' + domains_map_xor');
a_unif_disp = a_unif_disp(dispbox,dispbox);
figure; imagesc(a_unif_disp); colorbar; axis equal; caxis([-3 3]); 
title('logPrTAD - logPrBackground vs. Dixon TAD Mapping');
end

%a_ret = (domains_map).*a_diag + a_diag';
%a_ret = a_diag;
a_ret = a_gmm;
return;

%Used for rotated stuff. May use log here
a_display = a_diag;

%Calculte forward probability
%This works because all other are one's, and we have only positive diagonals
Pf = sum(log(a_prob),2); %row sum
Pb = sum(log(a_prob),1)'; %col sum

%plot(1:size(Pf),Pf,1:size(Pb),Pb);
%return

%Do Directionality index analisys
A = sum(a_diag,2);
B = sum(a_diag,1)';
E = (A+B)/2;
DI0 = ( (B-A)./abs(B-A) ) .* ( ((A-E).^2)./E + ((B-E).^2)./E );
DI1 = smooth(DI0,10); DI=DI1;%1
DI(DI>DI_TRIM)=DI_TRIM;
DI(DI<-DI_TRIM)=-DI_TRIM;
DI_SEQ = round(2*DI/DI_TRIM) + 3;
DI_SEQ(isnan(DI_SEQ)) = 3;

for i = 1:100
	likelystates = hmmviterbi(DI_SEQ, TRANS_GUESS, EMIS_GUESS);
	[TRANS_GUESS, EMIS_GUESS] = hmmestimate(DI_SEQ,likelystates);
end


%Prepare vectors for displaying
a_rot = imrotate(a_display,45);
a_width = size(a_rot,1);
diag_height = (MAX_DIAG-MIN_DIAG)/(2*sqrt(2));

di_plt = interp1(DI,linspace(1,numel(DI),a_width));
di_plt = (di_plt/max(di_plt(:)))*diag_height+(a_width/2+diag_height);

di_ref_x = [0 a_width];
di_ref_y = [a_width/2+diag_height a_width/2+diag_height];

di_likelystates = round(interp1(likelystates,linspace(1,numel(likelystates),a_width))-2)*diag_height + a_width/2 + diag_height;

%Plot forward/backword probability estimates
%normalize with the same maximum!
di_pf = interp1(Pf,linspace(1,numel(DI),a_width));
di_pf = smooth(di_pf,10);
di_pb = interp1(Pb,linspace(1,numel(DI),a_width));
di_pb = smooth(di_pb,10);

pf_norm = max(di_pf);

di_pf = di_pf/pf_norm*15*diag_height + (a_width/2-15*diag_height);
di_pb = di_pb/pf_norm*15*diag_height + (a_width/2-15*diag_height);

imagesc(a_rot);
colorbar;
axis equal;
hold on;
plot(di_plt,'Color','w','LineStyle','-');
plot(di_ref_x,di_ref_y,'Color','w','LineStyle',':');
plot(di_likelystates,'Color','r','LineStyle','-');
%plot(di_pf,'Color',[1 0 0],'LineStyle','-');
%plot(di_pb,'Color',[0 1 0],'LineStyle','-');

% 3 = down, 2 = no, 1 = up
% Looking for 3->1 transitions. Mark all 2's along the way
[startIdx,endIdx] = regexp(char(likelystates+48),'32*1');
a = zeros(size(likelystates));
for range = [startIdx; endIdx],
	a(range(1):range(2))=1;
end
plot((interp1(a,linspace(1,numel(likelystates),a_width))-2)*diag_height + a_width/2 - 0*diag_height,'Color','w','LineStyle','-');

hold off;
end
