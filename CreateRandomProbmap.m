%
% Maybe add noise
% Add another level of hierarchy
% Different strength of interaction
% Colorize using the interaction region only, not the entire box
%
function p = CreateRandomProbmap(sz)
	p = -ones(sz);

	%Basic TADs
	box = 1:(randi(35)+5);
	bounds = zeros(1,sz+1);
	bounds(1) = 0;
	while box(end) < sz
		bounds(box(end)+1) = 1;%Bound is start of next TAD
		p(box,box) = 1;
		box = (box(end)+1):(box(end)+randi(35)+5);
	end
	box = box(1):sz;
	bounds(end) = 0;
	p(box,box) = 1;

	%Level 2
	bnds = find(bounds);
	prevTAD = -1;
	for i = [1,bnds;bnds,sz+1];
		if prevTAD ~= -1
			fprintf('%d %d %d %d\r\n',prevTAD(1),prevTAD(2),i(1),i(2)-1);
			box = prevTAD(1):(i(2)-1);
			p(box,box) = p(box,box) + 0.25 + 0.25*rand(1);
			prevTAD = -1;
		elseif rand(1)>0.3
			prevTAD = i;
		end
	end
end
