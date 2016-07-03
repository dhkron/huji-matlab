% Can be used for other resolutions too
function [] = Liberman5kConvert()
	for i = 1:22
		Log(['Parsing chr' num2str(i)]);
		a = LibermanLoader(i,25000);
		Log();

		Log('Writing to disk');
		outPath = ['/cs/cbio/gil/raw/GSE63525_5K_Liberman/Processed/chr' num2str(i) '_25k.txt'];
		dlmwrite(outPath,full(a),'delimiter','\t');
		Log();
	end
