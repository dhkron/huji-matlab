%Assumes lower diagonal is all zeros or matrix is symmetric
function [c] = MakeSymmetric(a)
	if sum(sum(tril(a)))-trace(a) == 0
		b = triu(a).*(1-eye(size(a)));
		c = a + b';
	else
		c = a;
	end
end
