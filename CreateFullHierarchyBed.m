function = CreateFullHierarchyBed(prefix,res,fBedTemplate)
for chr = 1:19
	fprintf('Starting chr %d\r\n',chr);
	bedPath = sprintf(fBedTemplate,chr);
	CreateSinglehierarchyBed(prefix,res,chr,'',bedPath,'');
end

