% Input: cost 0 TADs
% Output: Hierarchy
function [] = TadTree(d,s,bnd,box_s,box_e,a_t,a_b)
	%Maximal merge distance, in blocks! assume this will be a paramete calculated with res and stuff
	maxMerge = 4000000/20000;

	%LLR height, in blocks! assume this magically got here. It is just 2mb/res
	H_llr = 2000000/20000;

	bounds = [0 bnd box_e-box_s+1];
	offset = box_s-1;

	textDraw = {};

	%fprintf('Range 1\t| Range2\t| V_DynProg\t| V_Super\n')
	%fprintf('-----------------------------------------------\n')
	import java.util.LinkedList;
	bl = LinkedList();

	%Generate the nice table of possible merges
	bound_matrix = zeros(numel(bounds)-2,8);
	bounds3 = [bounds(1:end-2) ; bounds(2:end-1) ; bounds(3:end) ];
	for i = 1:size(bounds3,2) 
		b = bounds3(:,i);
		s1 = b(1)+1;
		e1 = b(2);
		s2 = b(2)+1;
		e2 = b(3);
		d1 = d(s1,e1);
		d2 = d(s2,e2);
		sv = s(s1+offset,e2+offset);
		mdiff = d1+d2-sv;
		bound_matrix(i,:) = [s1 ,e1, s2, e2, d1,d2,sv,mdiff];
	
		xp = [s1 e1 e1 NaN];
		yp = [s1 s1 e1 NaN];
		patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
		textDraw{end+1} = [e1 s1 0];
	end
	if i>0 
		xp = [s2 e2 e2 NaN];
		yp = [s2 s2 e2 NaN];
		patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
		textDraw{end+1} = [e2 s2 0];
	end
	%bound_matrix

	mergeNumber = 1;
	while size(bound_matrix,1)>1
		sz = size(bound_matrix,1);
		min_diff = Inf;
		min_diff_s = -1;
		min_diff_e = -1;
		min_i = 0;
		%First, find the minimal merge
		for i = 1:sz
			tmerge_cost = bound_matrix(i,8);
			m_d1 = bound_matrix(i,2)-bound_matrix(i,1);
			m_d2 = bound_matrix(i,4)-bound_matrix(i,1);

			%Calculate amount of pixels in rectangle that are calculated with true LLR
			%LLR height = H
			%Left TAD base = D1
			%Right TAD base = D2
			%The diff is D1*D2 - partAboveH
			%partAboveH = big_triangle_above_h - the_two_smaller_triangles_above_h
			%Also we dont care about sqrt(2) factor
			S0 = m_d1*m_d2 ;
			bigTrBase = max(m_d1 + m_d2 - H_llr,0);
			smallTrBase1 = max(m_d1 - H_llr,0);
			smallTrBase2 = max(m_d2 - H_llr,0);
			S = S0 - (bigTrBase^2 - smallTrBase1^2 - smallTrBase2^2)/2; %/2 cause triangle
			tmerge_cost = tmerge_cost / S;
			%fprintf('%d-%d\t| %d-%d \t| %.1f\t| %.1f\t| %.1f\n',s1,e1,s2,e2,d_1+d_2,s_v,tmerge_cost);
			if tmerge_cost < min_diff && abs(bound_matrix(i,1)-bound_matrix(i,4)) < maxMerge && max(m_d1,m_d2) > 5
				min_diff = tmerge_cost;
				min_diff_s = bound_matrix(i,1);
				min_diff_e = bound_matrix(i,4);
				min_i = i;
			end
		end
		if min_diff == Inf
			break
		end
		
		%Update table entries
		%There are two possible steps, but sometimes only one is applicible
		if min_i > 1
			c = bound_matrix(min_i-1,:);
			
			s_p1 = min_diff_s + offset;
			s_p2 = min_diff_e + offset;
			s_p3 = c(1) + offset;
			new_s = s(s_p1, s_p2);
			big_s = s(s_p3, s_p2);
			if new_s == 0 || big_s == 0
				%fprintf('Has Zero! Fixing... \r\n');
				s(s_p1, s_p2) = ExtendSupersum(s_p1, s_p2, a_t, a_b);
				s(s_p3, s_p2) = ExtendSupersum(s_p3, s_p2, a_t, a_b);
				new_s = s(s_p1, s_p2);
				big_s = s(s_p3, s_p2);
			end

			new_mdf = c(5) + new_s - big_s;
			bound_matrix(min_i-1,:) = [c(1),c(2),c(3),min_diff_e,c(5),new_s,big_s,new_mdf];
		end
		if min_i < sz
			c = bound_matrix(min_i+1,:);

			s_p1 = min_diff_s + offset;
			s_p2 = min_diff_e + offset;
			s_p3 = c(4) + offset;
			new_s = s(s_p1, s_p2);
			big_s = s(s_p2, s_p3);
			if new_s == 0 || big_s == 0
				%fprintf('Has Zero! Fixing... \r\n');
				s(s_p1, s_p2) = ExtendSupersum(s_p1, s_p2, a_t, a_b);
				s(s_p2, s_p3) = ExtendSupersum(s_p2, s_p3, a_t, a_b);
				new_s = s(s_p1, s_p2);
				big_s = s(s_p2, s_p3);
			end

			new_mdf = new_s + c(6) - big_s;
			bound_matrix(min_i+1,:) = [min_diff_s,c(2),c(3),c(4),new_s,c(6),big_s,new_mdf];
		end

		xp = [min_diff_s min_diff_e min_diff_e NaN];%min_diff_e];
		yp = [min_diff_s min_diff_s min_diff_e NaN];%min_diff_s];
		patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);

		pos = (min_diff_s + min_diff_e)/2;
		
		textDraw{end+1} = [min_diff_e, min_diff_s, mergeNumber];	

		mergeNumber = mergeNumber + 1;

		bound_matrix = removerows(bound_matrix,'ind',min_i);
		%fprintf('%d-%d\n',min_diff_s,min_diff_e);
		
		%keyboard;
	end

	for c = textDraw
		drw = c{1};
		min_diff_e = drw(1);
		min_diff_s = drw(2);
		mergeNumber = drw(3);
		rad = 13;
		if mergeNumber == 0
			rad = 6;
			tsize = 1;
		elseif mergeNumber > 9
			tsize = 8;
		else
			tsize = 10;
		end
		rectangle('Position',[min_diff_e-rad/2 min_diff_s-rad/2 rad rad],'Curvature',[1,1,],'FaceColor','white');
		text('position',[min_diff_e, min_diff_s],'fontsize',tsize,'string',sprintf('%d',mergeNumber),'HorizontalAlignment','center','Color','black','FontName','arial');
	end
end
