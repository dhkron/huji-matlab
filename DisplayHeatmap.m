function DisplayHeatmap(a,range,box)
	figure('units','normalized','outerposition',[0 0 1 1])
	imagesc(a(box,box),range);
	axis equal;
	colorbar;
end
