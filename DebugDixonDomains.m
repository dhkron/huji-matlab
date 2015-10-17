function [] = DebugDixonDomains(fMatrix,fDomains,res)
	a = load(fMatrix);
	b = load(fDomains);
	DisplayHeatmap(log(a),0,0,'red');
	hold on;
	ax = axis;
	if ischar(res)
		res = str2num(res);
	end
	for c = b'
		c = c/res;
		i = c(1);
		j = c(2);
		plot([i,i],[ax(3),ax(4)],'g--');
		plot([j,j],[ax(3),ax(4)],'r-.');
	end
end
