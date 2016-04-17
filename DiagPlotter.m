% Plots diag-mean plots for all TAD-boundary boxes
% Produces nice graph that helps merging TADs
%
% Assumes
%	a = load('/tmp/test_em_matrix.txt');c = a(201:300,201:300);
% Test:
% d(12:73,12:73) = 1;
% d(x,y) = d(x,y)*(100 - abs(x-y) );
% DiagPlotter(d)

% Trick of tommy
% (D-BG)' \ (E-BG)' = regress((E-BG)',(D-BG)') 
function [boxes] = DiagPlotter(c,res)

bounds = [11,38,85,100,120,139,146]; %chr10 151:300 Gatto
bounds = [85,100,120,139,146]; %chr10 151:300 Gatto
bounds = [38,85,100,139,146]; %chr10 151:300 Gatto
bounds = [38,85,139,146]; %chr10 151:300 Gatto
bounds = [11,38,54,85,139]; %chr10 151:300 Gatto
bcp = [bounds(1:end-1) ; bounds(2:end) ];
figure;
colors = get(gca,'colororder');
legendInfo = {};
legendPoly = {};
patchs = {};
patchCode = 'A';%Actually 66 but nvm
i = 1;
boxes = {};
for b1idx = 1:size(bcp,2)
	b1 = bcp(:,b1idx);
	for b2idx = 1:size(bcp,2)
		b2 = bcp(:,b2idx);

		% patch, skip high order regions
		if abs(b1idx-b2idx)>1
			continue
		end

		% Different action if TAD or Overlap Region (or Hovering Overlap Region), regarding diags range
		x1 = b1(1);
		y1 = b2(1);
		x2 = b1(2)-1;
		y2 = b2(2)-1;
		% If points above diagonal
		if x1 <= y1 && x2 <= y2
			patchName = char(patchCode);

			t = c(x1:x2,y1:y2);
			[L,K] = DiagPlot(t);

			legendInfo{i} = sprintf('P1=%d,%d ; P2=%d,%d',x1,y1,x2,y2); % or whatever is appropriate
			% TAD region
			if x1 == y1 && x2 == y2
				L = L( (numel(L)+3)/2 : end );
				rng = 1:numel(L);
			% Overlap region
			elseif x2 == y1
				printf('Impossible!\r\n');
				L = L(2:end);
				rng = 1:numel(L);
			else
				rng = (y1-x2-1)+(1:numel(L));
			end
			L = smooth(L);
		

			[polystuff,errstuff] = polyfit(log(rng'*res),log(L),1);
			m = polystuff(1);
			logc=polystuff(2);
			R2 = 1 - errstuff.normr^2 / norm(log(L)-mean(log(L)))^2;
			legendPoly{i} = sprintf('%s: y=%g*x^{%g} R2=%g',patchName,exp(logc),m,R2);

			boxes{patchName} = {rng,L};

			pcolor = colors(1+mod(i-1,size(colors,1)),:);
			if x1 == y1 && x2 == y2
				patchs{i} = {[x1 x1 x2],[y1 y2 y2],pcolor,patchName};
			else
				patchs{i} = {[x1 x1 x2 x2],[y1 y2 y2 y1],pcolor,patchName};
			end
			patchCode = patchCode + 1;

			subplot(2,1,1);
			loglog(res*rng,L,'Color',pcolor);
			hold on;

			subplot(2,1,2);
			plot(res*rng,L);
			hold on;
			
			i = i + 1;
		end
	end
end
subplot(2,1,1);
legend(legendInfo,'Location','southwest');
grid on;
subplot(2,1,2);
legend(legendPoly,'Location','northeast');
grid on;

DisplayHeatmap(log2(c+1),0,0,'red');
hold on;
for p = patchs
	i = p{1};
	patch(i{2},i{1},i{3},'FaceAlpha',0.5,'EdgeColor','none');
	text(mean(i{2}),mean(i{1}),i{4},'FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center')
end

return
