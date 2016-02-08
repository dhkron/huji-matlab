function [a,a_llr,s] = FastCreateSingleHierarchyBed(prefix,res,chr,box,bedPath,figPath)
	
	hasBed = ( exist('bedPath','var') && numel(bedPath)>0 );	
	hasFig = ( exist('figPath','var') && numel(figPath)>0 );	

	template = '/cs/cbio/gil/htad-chain/output/%s.%s.chr%d.txt';
	a = (sprintf(template,prefix,'matrix',chr));
	a_llr = (sprintf(template,prefix,'llr',chr));
	s = (sprintf(template,prefix,'supersum',chr));

	if ~hasBed
		bedPath = '';
	end
	if ~hasFig
		figPath = '';
	end

	%Passing chrnumber is stinn nesessary. need to remove the 'chr' part before passig it to matlab
	CreateSingleHierarchyBed(a,a_llr,s,prefix,res,chr,box,bedPath,figPath);
end
