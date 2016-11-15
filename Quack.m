%Height is 5M - @TadTree
function [a,b,c] = Quack(chrNum,prefix)

chrNum = 17;
prefix = 'mCO';
matPath = sprintf('/cs/cbio/gil/htad-chain/output/%s.matrix.chr%s.txt',prefix,chrNum);
bedPath = sprintf('/cs/cbio/gil/htad-chain/output/bed/chr%s_%s_full.bed',chrNum,prefix);
dxnPath = sprintf('/cs/cbio/gil/raw/Dixon/%s_Domains/combined/chr%s.combined.domain',prefix,chrNum);

bedOutPath = '/tmp/deleteme';
[a,b,c] = RebuildMatrices(matPath,bedPath,dxnPath, bedOutPath);

end
