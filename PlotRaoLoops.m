function PlotRaoLoops(H,x1,x2,y1,y2,offset,res)
	x1 = x1/1000+1;
	x2 = x2/1000+1;
	y1 = y1/1000+1;
	y2 = y2/1000+1;
	box1=x1-offset:x2+offset;
	box2=y1-offset:y2+offset;
	%If you are not sure what H is, look in README located in raw Rao/Liberman data
	DisplayHeatmap(log2(1+H(box1,box2)),[],[],'red')
	%Pos1-box1 is effectivly offset.
	w = x2-x1;
	h = y2-y1;
	%Black offset magic
	rectangle('Position',[offset+1, offset+1, w, h])
end
