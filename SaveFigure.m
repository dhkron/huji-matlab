function SaveFigure(handle, path)
	handle.Units = 'centimeters';
	handle.PaperUnits = 'centimeters';
	handle.PaperPosition = handle.Position;
	saveas(handle,path); %In matlab, the extension matters
end
