function [tad_size, tad_color] = CalcTADStrength(chrNumber,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE,doPlot)

domains = load(sprintf('TADs/domains.IMR90.chr%d', chrNumber));
domains = floor(domains/BLOCK_SIZE)+1; %To start from 1
tad_size = zeros(1,a_size);
tad_color = zeros(1,a_size);
tad_index = 1
for i = MIN_DIAG:MAX_DIAG
	for domain = domains'
		domain(2) = min(domain(2),a_size); %Fix overflow

		tad_sum = sum(sum((a(domain(1):domain(2),domain(1):domain(2))))); %Can use log here
		tad_count = numel(a(domain(1):domain(2),domain(1):domain(2)));
		tad_size(tad_index) = domain(2)-domain(1);
		tad_color(tad_index) = tad_sum/tad_count;
		tad_index = tad_index + 1;
	end
end

%Plot TAD strength vs. Length
if exist('doPlot','var') && doPlot
	figure;
	scatter(tad_size*BLOCK_SIZE,tad_color);
	title('Mean interaction strength vs. TAD length');
	xlabel('Length (bp)');
	ylabel('Mean Interaction strength');
	set(gca,'xscale','log');
	set(gca,'yscale','log');
	coeffs = polyfit(log(tad_size*BLOCK_SIZE),log(tad_color),1);
	hold on;
	plot(tad_size*BLOCK_SIZE,((tad_size*BLOCK_SIZE).^coeffs(1))*exp(coeffs(2)));
	hold off;
	legend('Original Data',sprintf('Line Fit ae^b\r\na=%d\r\nb=%d',exp(coeffs(2)),coeffs(1)));
end

end
