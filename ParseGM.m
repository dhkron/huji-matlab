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
	sigmaBg(i_cell) = sqrt(db.Sigma(i_bg));
	sigmaIn(i_cell) = sqrt(db.Sigma(i_in));
	probBg(i_cell) = db.ComponentProportion(i_bg);
	probIn(i_cell) = db.ComponentProportion(i_in);
end % Loop

if exist('plotLog','var') && plotLog
	x_vals = (1:numel(a)) * 40000;
	a_bg = polyfit(log(x_vals), log(meansBg),1);
	a_in = polyfit(log(x_vals), log(meansIn),1);
	%fprintf('Power law of background: %f (%f)\r\n',a_bg)
	%fprintf('Power law of interaction: %f (%f)\r\n',a_in);

	figure;
	hold on;
	loglog(x_vals, meansBg,'bx-','LineWidth',2);
	hold on;
	loglog(x_vals, meansIn,'rx-','LineWidth',2);
	loglog(x_vals, meansBg+sigmaBg,'b--');
	loglog(x_vals, meansBg-sigmaBg,'b--');
	loglog(x_vals, meansIn+sigmaIn,'r--');
	loglog(x_vals, meansIn-sigmaIn,'r--');
	%plot(x_vals,x_vals.^a_bg(2)*exp(a_bg(2)),'b-');
	%plot(x_vals,x_vals.^a_in(1)*exp(a_in(2)),'r-');
	legend('Background','Interaction')
	title('Mean (with sigma) vs. Block distance')
	xlabel('Base pairs distance')
	ylabel('Mean +/- sigma')

end

end % Function

