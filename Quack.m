%Height is 5M - @TadTree
function [a,b,c] = Quack(chrNum,prefix)

chrNum = 17;
prefix = 'mCO';
matPath = sprintf('/cs/cbio/gil/htad-chain/output/%s.matrix.chr%d.txt',prefix,chrNum);
bedPath = sprintf('/cs/cbio/gil/htad-chain/output/bed/chr%d_%s_full.bed',chrNum,prefix);
dxnPath = sprintf('/cs/cbio/gil/raw/Dixon/%s_Domains/combined/chr%d.combined.domain',prefix,chrNum);

bedOutPath = '/tmp/deleteme';
[a,b,c] = RebuildMatrices(matPath,bedPath,dxnPath, bedOutPath);

end
