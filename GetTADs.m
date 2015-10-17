function [a_b,a_t,a_c] = GetTADs(chrNumber,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE)

domains = load(sprintf('TADs/domains.IMR90.chr%d', chrNumber));
[a_b,a_t,a_c] = GetTADsMask(domains,MIN_DIAG,MAX_DIAG,a_size,BLOCK_SIZE);

end
