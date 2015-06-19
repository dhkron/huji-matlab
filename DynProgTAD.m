function [dynmap, dynsol] = DynProgTAD(supermapraw,s,e, pyrraw)
supermap = supermapraw(s:e,s:e);
pyr = pyrraw(s:e,s:e);
dynmap = zeros(size(supermap));
dynsol = zeros(size(supermap));

supermapfl = fliplr(supermap);

%perms = {};
%for i=1:(e-s+1)
%	perms{i} = randperm(i);
%end

a_size = size(supermap,1);
for dg = 0:a_size
	for pos = 1:(a_size-dg)
		i = pos;%row
		j = dg+pos;%col

		if abs(i-j)<5
			dynmap(i,j) = supermap(i,j);
		else
			ptad = supermap(i,j);
			midtype = 0;
			%rng = i:j-1;
			%rng = rng(perms{numel(rng)});
			%for mid = rng
			for mid = i:j-1
				flipj = a_size+1-j;
				v = dynmap(i,mid) + dynmap(mid+1,j);
				%v = v + 10*(j-i)*( - pyr(i,j) + pyr(i,mid) + pyr(mid+1,j) );
				% + Interaction square value?
				% How does -100 affect the structure?
				%If I add a measure of the cross-tad interaction?
				%How can it find DISTANT interacting TADs?
				if v > ptad
					ptad = v;
					midtype = mid;
				end
			end
			dynmap(i,j) = ptad;
			dynsol(i,j) = midtype;
		end
	end
end
%figure; imagesc(dynmap); axis equal; colorbar;

disp('Finished! Printing solution...')
printSolution(dynsol,1,e-s+1,0); %offset=s-1
fprintf('\r\n');

end

function printSolution(dynsol,i,j,offset)
	current = dynsol(i,j);
	if current == 0
		fprintf('%d-%d',i+offset,j+offset)
		if abs(i-j) > 5
			hold on;
			ax = axis;
			plot([i,i],[ax(3),ax(4)],'b--');
			plot([ax(1),ax(2)],[i,i],'b--');
			plot([j,j],[ax(3),ax(4)],'r-.');
			plot([ax(1),ax(2)],[j,j],'r-.');
			xp = [i i j j];
			yp = [i j j i];
			patch(xp,yp,'magenta','FaceAlpha',0.25,'EdgeColor','none');
		end
	else
		fprintf('( ')
		printSolution(dynsol,i,current,offset);
		fprintf(' , ')
		printSolution(dynsol,current+1,j,offset);
		fprintf(' )')
	end
end
