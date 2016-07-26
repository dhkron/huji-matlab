% Can be used for other resolutions too
function [] = Liberman5kConvert(prefix, path)
	for i = 1:20
		Log(['Parsing chr' num2str(i)]);
		a = LibermanLoader(i,25000,path);
		Log();

		Log('Writing to disk');
		outPath = ['/cs/cbio/gil/raw/GSE63525_5K_Liberman/Processed/' prefix '.chr' num2str(i) '_25k.txt'];
		dlmwrite(outPath,MakeSymmetric(full(a)),'delimiter','\t');
		Log();
	end
