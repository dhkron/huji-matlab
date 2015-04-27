function [a_ret] = TAD(chrNumber)
DATA_FILE = sprintf('HiC-CSV-Matrices/IMR90_chr%d_40k', chrNumber);
MIN_DIAG = 1;%0 is the middle one
MAX_DIAG = 55;%55
BLOCK_SIZE = 40000;

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

%Quack
a_ret = 0;

%Load
fprintf('Loading... ');
a=load(DATA_FILE);
fprintf('Done\r\n');

%Normalize
if 1
	fprintf('Normalizing... ');
	a = a + ones(size(a));
	for i=1:100
		a = (a./repmat(sum(a),size(a,1),1))';
	end
	fprintf('Done\r\n')
end

%Trim diagonal
diag_matrix = zeros(size(a));
for x = MIN_DIAG:MAX_DIAG
	diag_matrix = diag_matrix + diag(ones(size(a,1)-x,1),x);
end
clear_diag = logical(diag_matrix);
a_diag = a .* clear_diag;

%Generate GMM data
a_size = size(a,1);
a_ret = {};
fprintf('Diagonalizing... ');
warning('off','stats:gmdistribution:FailedToConverge');
for i = MIN_DIAG:MAX_DIAG
	dgn = diag(a_diag,i);
	dgn_clean = dgn(dgn~=0);

	GMModel = fitgmdist(dgn_clean,2);
	a_ret{i}=GMModel;
end
w = warning('query','last');
fprintf('Done\r\n');

%Use public TAD data to extract intraTAD stuff
fprintf('Extracting TADs... ');
domains = load(sprintf('TADs/domains.IMR90.chr%d', chrNumber));
domains = floor(domains/BLOCK_SIZE)+1;
domains_map = zeros(a_size);
for i = MIN_DIAG:MAX_DIAG
	start_pos = a_size*i + 1;
	domains_diag = zeros(1,a_size-i);
	for domain = domains'
		bnd_start = domain(1)+1;
		bnd_end = domain(2)-i;
		min_bnd = min(bnd_start,bnd_end);
		min_bnd = max(1, min_bnd);
		max_bnd = max(bnd_start,bnd_end);
		max_bnd = min(max_bnd, a_size-i);
		domains_diag(min_bnd:max_bnd) = ~domains_diag(min_bnd:max_bnd);

		%min_bnd = max(1, bnd_start);
		%max_bnd = min(bnd_end, a_size-i);
		%domains_diag(min_bnd:max_bnd) = 1; %Do not extend TADs to tiles
	end
	domains_map(start_pos:a_size+1:end) = domains_diag;
	%domains_map(start_pos:a_size+1:end) = round(rand(1,a_size-i)); %Random model
end
a_tad = domains_map .* a_diag;
a_bound = (1-domains_map) .* a_diag; %For bounderies
fprintf('Done\r\n');

%Extract probability for TAD diagonals only
figure;
iplt = 1;
for i = [2,5,10,20,30,45]
	dgn = log(1+diag(a_tad,i)); %This does not do anything apparently
	%dgn = diag(a_tad,i);
	dgn_clean = dgn(dgn~=0);
	[f,xi] = ksdensity(dgn_clean,dgn,'function','pdf');

	subplot(2,3,iplt);
	scatter(xi,f,11,'filled');
	hold on;

	dgn = log(1+diag(a_bound,i)); 
	%dgn = diag(a_bound,i);
	dgn_clean = dgn(dgn~=0);
	[f,xi] = ksdensity(dgn_clean,dgn,'function','pdf');

	scatter(xi,f,11,'filled');
	title(sprintf('Diagonal=%d',i));
	legend('intra-domains','bounderies');
	iplt = iplt+1;
end
suptitle('Probability density of interactions strength, for several diagonals');
a_ret = a_tad;
%figure; imagesc(domains_map); axis equal; colorbar;
return;

%Now, remember how TADs are looking in the 2d unrotated map.
%I should add 1 to the START of the TAD only! Depends on stuff. Just add to begginig and see how it works.
%When they collide - reverse? or just always put lower and higher? This might word better
%Then, the world!

%Set return value to be diagonalized matrix
a_ret = a_diag;
%return

%Generate logscale (disabled)
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
