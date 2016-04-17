% This function takes the structure produced from DiagPlotter.m
% Takes Background model and compares to Tad-Tad interaction
% 
% It parses automatically the boxes data structure

function DiagPlotterHrrcAnalyzer(boxes,a,b,bg)
	A = boxes{a}{2};
	B = boxes{b}{2};
	BG = bg';
	l = min([numel(A), numel(B), numel(bg)])
	regress(B(1:l)-BG(1:l), A(1:l)-BG(1:l))
