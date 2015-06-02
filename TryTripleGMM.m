function TryTripleGMM(a_diag,MIN_DIAG,MAX_DIAG)

%Generate GMM data
fprintf('GMMing (3 Param)... ');
a_gmm = AnalyzeGMM(a_diag,MIN_DIAG,MAX_DIAG,3);
fprintf('Done\r\n');

figure;
for i = MIN_DIAG:MAX_DIAG
	x = i*BLOCK_SIZE;
	loglog(x,a_gmm{i}.mu(1),'*r');
	hold on;%for loglog
	loglog(x,a_gmm{i}.mu(2),'*g');
	loglog(x,a_gmm{i}.mu(3),'*b');
	title('Mean'); xlabel('Distance (bp)');
end
figure;
for i = MIN_DIAG:MAX_DIAG
	x = i*BLOCK_SIZE;
	loglog(x,a_gmm{i}.Sigma(1),'*r');
	hold on;%for loglog
	loglog(x,a_gmm{i}.Sigma(2),'*g');
	loglog(x,a_gmm{i}.Sigma(3),'*b');
	title('Sigma'); xlabel('Distance (bp)');
end
figure; hold on;
for i = MIN_DIAG:MAX_DIAG
	x = i*BLOCK_SIZE;
	plot(x,a_gmm{i}.ComponentProportion(1),'*r');
	plot(x,a_gmm{i}.ComponentProportion(2),'*g');
	plot(x,a_gmm{i}.ComponentProportion(3),'*b');
	title('Posterior'); xlabel('Distance (bp)');
end
end
