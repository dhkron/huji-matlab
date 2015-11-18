function SaveFigure(handle, path)
	handle.Units = 'centimeters';
	handle.PaperUnits = 'centimeters';
	handle.PaperPosition = handle.Position;

	[~,fname,fext] = fileparts(path);

	Log(['Saving file "',fname,fext,'"']);
	saveas(handle,path); %In matlab, the extension matters
	Log();
end
