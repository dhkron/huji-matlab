function [f] = DebugDixonDomains(fMatrix,fDomains,res,box)
	if ischar(fMatrix) 
		a = load(fMatrix);
	else
		a = fMatrix;
	end
	if ischar(fDomains)
		b = load(fDomains);
	else
		b = fDomains;
	end
	hasBox = ( exist('box','var') && numel(box)>1);
	if hasBox
		f = DisplayHeatmap(log(a),0,box,'orange');
	else
		f = DisplayHeatmap(log(a),0,0,'orange');
	end
	hold on;
	ax = axis;
	if ischar(res)
		res = str2num(res);
	end
	for c = b'
		c = c/res;
		i = c(1);
		j = c(2);
		if i > box(1) && i < box(end)
			i = i-box(1)+1;
			plot([i,i],[ax(3),ax(4)],'g--');
			plot([ax(1),ax(2)],[i,i],'g--');
		end
		if j > box(1) && j < box(end)
			j = j-box(1)+1;
			plot([j,j],[ax(3),ax(4)],'b-.');
			plot([ax(1),ax(2)],[j,j],'b-.');
		end
	end
end
