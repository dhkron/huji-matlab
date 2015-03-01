function [] = TAD(chrNumber)
DATA_FILE = sprintf('HiC-CSV-Matrices/IMR90_chr%d_40k', chrNumber);
MIN_DIAG = 5;
MAX_DIAG = 55;

DI_TRIM = 1.0/250;
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
fprintf('Normalizing... ');
for i=1:100
	a = (a./repmat(sum(a),size(a,1),1))';
	a(isnan(a))=0;
end
fprintf('Done\r\n')

%Probability vector
%Should derive here probability model.
%However, i've lost it when I normalized.
%Should derive probabilistic model from diag histogram,
%And calculate probabilities using it.
%It is ok that PDF is larger than 1,
%It seems that I need the CDF?
%PDF makes no sense in continiues?
diag_model = zeros(1,MAX_DIAG);
a_size = size(a,1);
for dn = MIN_DIAG:MAX_DIAG
	diag_model(dn) = sum(diag(a,dn-1));
end
figure;
hold on;
for i = MIN_DIAG:MAX_DIAG
	dgn = diag(a,i);
	dgn = dgn(dgn~=0);
	[f,xi] = ksdensity(dgn,'function','pdf');
	plot(xi,f);
end
hold off;
%figure;plot(diag_model);
return

%Trim diagonal
diag_matrix = zeros(size(a));
for x = MIN_DIAG:MAX_DIAG
	diag_matrix = diag_matrix + diag(ones(size(a,1)-x,1),x);
	%diag_matrix = diag_matrix + diag(ones(size(a,1)-x,1),-x);
end
clear_diag = logical(diag_matrix);
a_diag = a .* clear_diag;

%Generate logscale (disabled)
a_log = a;

%Do Directionality index analisys
A = sum(a_diag,2);
B = sum(a_diag,1)';
E = (A+B)/2;
DI0 = ( (B-A)./abs(B-A) ) .* ( ((A-E).^2)./E + ((B-E).^2)./E );
DI1 = smooth(DI0,10); DI=DI0;%1
DI(DI>DI_TRIM)=DI_TRIM;
DI(DI<-DI_TRIM)=-DI_TRIM;
DI_SEQ = round(2*DI/DI_TRIM) + 3;
DI_SEQ(isnan(DI_SEQ)) = 3;

for i = 1:100
	likelystates = hmmviterbi(DI_SEQ, TRANS_GUESS, EMIS_GUESS);
	[TRANS_GUESS, EMIS_GUESS] = hmmestimate(DI_SEQ,likelystates);
end


%Prepare vectors for displaying
a_rot = imrotate(a_log,45);
a_width = size(a_rot,1);
diag_height = (MAX_DIAG-MIN_DIAG)/(2*sqrt(2));

di_plt = interp1(DI,linspace(1,numel(DI),a_width));
di_plt = (di_plt/max(di_plt(:)))*diag_height+(a_width/2+diag_height);

di_ref_x = [0 a_width];
di_ref_y = [a_width/2+diag_height a_width/2+diag_height];

di_likelystates = round(interp1(likelystates,linspace(1,numel(likelystates),a_width))-2)*diag_height + a_width/2 + diag_height;

imagesc(a_rot,[0 0.01]);
colorbar;
hold on;
plot(di_plt,'Color','w','LineStyle','-');
plot(di_ref_x,di_ref_y,'Color','w','LineStyle',':');
plot(di_likelystates,'Color','r','LineStyle','-');

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
