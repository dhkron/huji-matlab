% Input: .bed file, heatmap, resolution
% Output: Hierarchy
function [B] = DrawBed(res,chrNum,prefix)
	%Bed is currently ignored since all beds are well formatted (yey for good coding (this is not good coding))
	fBed = fopen(sprintf('/cs/cbio/gil/htad-chain/output/bed/chr%d_%s_full.bed',chrNum,prefix),'r');
	[B,C] = textscan(fBed,'chr%d\t%d\t%d\t%s\t%f',Inf);
	fclose(fBed);

	a = load(sprintf('/cs/cbio/gil/htad-chain/output/%s.matrix.chr%d.txt',prefix,chrNum));
	box = 1:size(a,2);
	DisplayHeatmap(log2(a+1),[0,6],box,'red');
	title(sprintf('Hierarchy TAD numbered nesting - %s chr%d',prefix,chrNum));
	hold on;

	BLen = size(B{1},1);
	for i = 1:BLen
		s = B{2}(i)/res;
		e = B{3}(i)/res;
		xp = [s,e,e,NaN];
		yp = [s,s,e,NaN];
		patch(xp,yp,'magneta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
	end
	for i = 1:BLen
		s = B{2}(i)/res;
		e = B{3}(i)/res;
		xp = [s,e,e,NaN];
		yp = [s,s,e,NaN];
		mergeNumber = i;
		rad = 13;
		if mergeNumber == 0
			rad = 6;
			tsize = 1;
		elseif mergeNumber > 9
			tsize = 8;
		else
			tsize = 10;
		end
		%Yes, in matlab you draw circle with rectangle
		rectangle('Position',[e-rad/2 s-rad/2 rad rad],'Curvature',[1,1,],'FaceColor','white');
		%text('position',[min_diff_e, min_diff_s],'fontsize',tsize,'string',sprintf('%d',mergeNumber),'HorizontalAlignment','center','Color','black','FontName','arial');
	end
	path = sprintf('/cs/grad/gilr02/Dropbox/Bio/NextMeeting/Hierarchies_%s_chr%d_Full.png',prefix,chrNum);
	g = gcf;
	SaveFigure(g,path);
	close(g);
end
