% Was made to look at TADs with bilinear line
% My theory was that bilinears are simply two tads that were mis-marked
% It looked to me like itwas the case, but Tommy was not convinced
function [] = MarginDisplayBi(c,x,y)
	box=floor(x/4):floor(y/4);
	box2=[box(1)-20:box(1)-1 box box(end)+1:box(end)+20];
	boxw = box2(end)-box2(1)+1;
	DisplayHeatmap(log2(1+c),[],box2,'red');
	rectangle('Position',[20 20 boxw-40 boxw-40]);
