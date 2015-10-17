function f = DisplayHeatmap(a,range,box,scheme)
	f = figure('units','normalized','outerposition',[0 0 1 1]);
	hasRange = ( exist('range','var') && numel(range)>1 );
	hasBox = ( exist('box','var') && numel(box)>1);
	hasScheme = exist('scheme','var') && numel(scheme)>0;
	if hasRange && hasBox
		imagesc(a(box,box),range);
	elseif hasRange
		imagesc(a,range);
	elseif hasBox
		imagesc(a(box,box));
	else
		imagesc(a);
	end
	if hasScheme
		red = [ones(256,1),(256:-1:1)'/256,(256:-1:1)'/256];
		colormap(red);
	end
	axis equal;
	colorbar;
end
