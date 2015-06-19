function DisplayHeatmap(a,range,box)
	figure('units','normalized','outerposition',[0 0 1 1])
	if exist('range','var')
		if exist('box','var')
			imagesc(a(box,box),range);
		else
			imagesc(a,range);
		end
	else
		imagesc(a);
	end
	axis equal;
	colorbar;
end
