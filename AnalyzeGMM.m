function a_gmm = AnalyzeGMM(a,MIN_DIAG,MAX_DIAG)
	a_gmm = {}
	warning('off','stats:gmdistribution:FailedToConverge');
	for i = MIN_DIAG:MAX_DIAG
		dgn = diag(a,i);
		dgn_clean = dgn(dgn~=0); %Zeros are bad for fitting?
		dgn_clean = log(dgn_clean+1); %Log-normal %Should replace 1 with something smarter

		options = statset('MaxIter',1000);
		GMModel = fitgmdist(dgn_clean,2,'Options',options);
		a_gmm{i}=GMModel;

		[~,msgid] = lastwarn;
		if strcmp(msgid, 'stats:gmdistribution:FailedToConverge')
			fprintf('Failed to converge at diag %d\r\n',i);
			warning('')
		end
	end
end
