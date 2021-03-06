%
% I guess pii should be 1000 or so
%
function [] = CreateAndSaveLotsOfImages(res,pixels_in_image) 
	path = '/cs/cbio/gil/htad-chain/output/';
	a_tmp = strcat(path,'m%s.matrix.chr%s.txt');
	s_tmp = strcat(path,'m%s.supersum.chr%s.txt');
	llr_tmp = strcat(path,'m%s.llr.chr%s.txt');
	d_tmp = '/cs/cbio/gil/raw/Dixon/m%s_Domains/combined/chr%s.combined.domain';

	%chrsz = load('/cs/cbio/gil/raw/ChrSize/mm9.txt');
	%chrsz = fliplr(chrsz); %So that chrsz(i) = chrsz(i,1) = {Chromsome-i}'s size [also chrsz(i,2) become chr number actually]

	for i = [1:19]
		for tp = {'ES','CO'}
			fprintf('Loading %s chr%s\r\n',tp{1},i);
			a_path = sprintf(a_tmp,tp{1},i);
			s_path = sprintf(s_tmp,tp{1},i);
			llr_path = sprintf(llr_tmp,tp{1},i);
			d_path = sprintf(d_tmp,tp{1},i);
			for checkf = {a_path,s_path,llr_path,d_path}
				if ~exist(checkf{1},'file')
					fprintf('%s does not exist\r\n',checkf{1});
				end
			end
			a = load(a_path);
			s = load(s_path);
			llr = load(llr_path);
			d = load(d_path);
			
			chrsz_i = size(a,1);%In blocks, of current, used to be chrsz(i)
			%BEFORE USING THIS CODE CHANGE DISPLAYHEATMAP TO MAKE INVISIBLE UNICORNS. Not required in "-nodisplay"
			fprintf('Saving %s chr%s\r\n',tp{1},i);
			box = 1:min(pixels_in_image,chrsz_i);
			box_index = 1;
			box_inc = floor(pixels_in_image/2);

			while numel(box) > 0
				if numel(box) < pixels_in_image && box_index > 1
					box = (box(end)-pixels_in_image+1):box(end);
					flagBreak = 1;
				else
					flagBreak = 0;
				end

				%Save figure of my TAD mapping
				g = DisplayHeatmap(log2(a+1),[0 6],box,'red');
				[d1 d2 d3] = DynProgTAD(s,box(1),box(end),0);

				%title(sprintf('My Method, %s Chr%d, Range:[%d-%d]bp',tp{1},i,box(1),box(end)));
				%SaveFigure(g,sprintf('/cs/cbio/gil/htad-chain/output/png/m%s.chr%s.part%d.mymethod.png',tp{1},i,box_index));

				%Add hierarchy layer and save
				TadTree(d1,s,find(d3),box(1),box(end),llr,-llr,res,i);
				title(sprintf('Hierarchy TAD numbered nesting - %s chr%s [blocks %d-%d]',tp{1},i,box(1),box(end)));
			
				%Update labels. LABELS MUST BE UPDATED AFTER ***LAST*** CHANGE TO GRAPH	
				newlbl = get(gca,'XTickLabel');
				for q = 1:numel(newlbl)
					newlbl{q} = num2str((str2num(newlbl{q})+box(1)-1)*(res/1e6));
				end
				set(gca,'XTickLabel',newlbl);
				xlabel('in megabase');
				newlbl = get(gca,'YTickLabel');
				for q = 1:numel(newlbl)
					newlbl{q} = num2str((str2num(newlbl{q})+box(1)-1)*(res/1e6));
				end
				set(gca,'YTickLabel',newlbl);
				ylabel('in megabase');

				SaveFigure(g,sprintf('/cs/cbio/gil/htad-chain/output/png/m%s.chr%s.part%d.hierarchy.png',tp{1},i,box_index));
				close(g);

				%Save figure of dixon results
				%g = DebugDixonDomains(a,d,40001,box);
				%title(sprintf('Dixon Published Results, %s Chr%d, Range:[%d-%d]bp',tp{1},i,box(1),box(end)));
				%SaveFigure(g,sprintf('/cs/cbio/gil/htad-chain/output/png/m%s.chr%s.part%d.dixonpub.png',tp{1},i,box_index));
				%close(g);

				box_index = box_index + 1;
				box = (box(1)+box_inc) : min(box(end)+box_inc,chrsz_i);

				if flagBreak == 1
					break
				end
			end
		end
	end
end
