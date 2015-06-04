function DisplayHeatmap(a,range,box)
	figure;
	imagesc(a(box,box),range);
	axis equal;
	colorbar;
end
