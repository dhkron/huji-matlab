function [] = BedFindPowerlaw(pathToBed,a,res)
T = readtable('/cs/cbio/gil/htad-chain/output/bed/matlab_compatible/chr14.txt','Delimiter',' ');
a = load('/cs/cbio/gil/htad-chain/output/mCO.matrix.chr14.txt');
res = 40000;
N = 250;
Skip = 3;

tads = containers.Map([0],[0]);

if 0
DIST = (Skip+1):N;
else 
base = 1.4;
DIST = unique(floor([0 base.^[0:log(N)/log(base)] N]));
end

XX = res*DIST;
YY = ones(N,1);

%1 if using exponential
flg = (DIST(end)-DIST(end-1)~=1);

for ind = 1:length(T.Type)
	%fprintf('%d %d\r\n',T.Start(ind),T.End(ind));
	if strcmp(T.Type(ind),'M')
		box1 = tads(T.End(ind))/res:T.End(ind)/res; %B2 to C = boxX, also RIGHT TAD
		box2 = T.Start(ind)/res:tads(T.Start(ind))/res; %A to B1 = boxY, also LEFT TAD
	
		box = (T.Start(ind)/res):(T.End(ind)/res);
		%fprintf('%d-%d %d-%d *** ',box2(1),box2(end),box1(1),box1(end));

		tads(T.End(ind)) = T.Start(ind);
		tads(T.Start(ind)) = T.End(ind);
		
		%prompt = '0 = nothing, 1 = print, 2 = quit ? ';
		xinp = 1;%input(prompt);
		
		if xinp == 1
			abox = a(box2,box1);
			offset = 1-size(abox,1);

			%DisplayHeatmap(log2(a+1),0,box(1)-10:box(end)+10,'red');caxis([0 5]);	
			%DisplayHeatmap(log2(abox+1),0,0,'red');caxis([0 5]);axis equal;title(sprintf('Merge area %d',ind));

			%fprintf('~~ %d-%d %d-%d *** ',box2(1),box2(end),box1(1),box1(end));
			%fprintf('Break at image %d *** offset %d ~~\r\n',box2(end)-box2(1)+10,offset);
		elseif xinp == 2
			break
		else 
			continue
		end
	else
		tads(T.Start(ind)) = T.End(ind);
		tads(T.End(ind)) = T.Start(ind);

		box = (T.Start(ind)/res):(T.End(ind)/res);
		box1 = box;
		box2 = box;
		abox = a(box,box);

		offset = 0;
	end

	[I,J] = find(~isnan(abox));
	L = NaN*ones(1,length(DIST));
	for i = (1+flg):length(DIST);
		%Tad mode: J-I==0 is the base diagonal
		%Merge mode: J-I==(1-size(abox,1)) is base diagonal (max is 1-size2)
		%That means, J-I==DIST(1) should translate to J-I==DIST(1) + (1-size1)
		if flg
			D = abox(J-I>DIST(i-1)+offset & J-I<=DIST(i)+offset);
		else
			D = abox(J-I==DIST(i)+offset);
		end
		if length(D) < 9, continue; end
		L(i) = nanmean(D);
	end

	II = find(~isnan(L));
	[b,bin,~,~,stat]=regress(log(L(II))', [log(XX(II)'), YY(II)]);
	R2 = stat(1);

	YH=exp( [log(XX(II)'), YY(II)] * b );
	fprintf('%s\t%d\t%.2f\t%.2f\t%.2f\r\n',T.Type{ind},length(box),b(1),b(2),R2);
	
	if T.Type{ind} == 'M' %isinf(b(1))
		%figure; loglog(XX(II),L(II),'o',XX(II),YH,'r-',XX(1:Skip),L(1:Skip),'go'); title(sprintf('%s%d',T.Type{ind},ind));
	end
	if b(1) > -0.7 
		%DisplayHeatmap(log2(a+1),0,box(1)-10:box(end)+10,'red');
	end
end

end
