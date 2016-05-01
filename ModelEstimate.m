function [] = ModelEstimate(matPath,bedPath,bedPathOut,emMatOut)

	[RES_Hrr,~,~] = RebuildMatrix(matPath,bedPath,0,bedPathOut);

	Log('Writing model estiamted matrix');
	dlmwrite(emMatOut,RES_Hrr);
	Log();

end
