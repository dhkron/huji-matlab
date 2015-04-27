function [meansBg, meansIn, sigmaBg, sigmaIn, probBg, probIn] = ParseGM(a, plotLog)
meansBg = zeros(size(a));
meansIn = zeros(size(a));
sigmaBg = zeros(size(a));
sigmaIn = zeros(size(a));
probBg = zeros(size(a));
probIn = zeros(size(a));
for i_cell = 1:numel(a)
	db = a{i_cell};
	if db.mu(1) < db.mu(2) %Higher pixel values are interactions
		i_bg = 1;
		i_in = 2;
	else
		i_bg = 2;
		i_in = 1;
	end
	meansBg(i_cell) = db.mu(i_bg);
	meansIn(i_cell) = db.mu(i_in);
	sigmaBg(i_cell) = db.Sigma(i_bg);
	sigmaIn(i_cell) = db.Sigma(i_in);
	probBg(i_cell) = db.ComponentProportion(i_bg);
	probIn(i_cell) = db.ComponentProportion(i_in);
end % Loop

a_bg = polyfit(log(1:numel(a)), log(meansBg),1);
a_in = polyfit(log(1:numel(a)), log(meansIn),1);
fprintf('Power law of background: %f (%f)\r\n',a_bg)
fprintf('Power law of interaction: %f (%f)\r\n',a_in);

if exist('plotLog','var') && plotLog
	loglog(1:numel(a), meansBg);
	hold on;
	loglog(1:numel(a), meansIn);
end
end % Function

