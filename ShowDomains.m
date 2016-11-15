function [domains] = ShowDomains(chrNumber)
	hold on;
	domains = load(sprintf('TADs/domains.IMR90.chr%s',chrNumber));
	for a = domains'
		b = a/40000*sqrt(2);
		line([b(1) b(1)],[-1e6 1e6],'Color',[0 1 0 0.5],'LineWidth',2);
		line([b(2) b(2)],[-1e6 1e6],'Color',[1 0 0 0.5],'LineWidth',1);
	end
	hold off;
end
