function CloneCNA_plot_normalized_results(Datafile,resultsfile,plotsdir,barcode)
%this function is used to plot 22 figures given the data, all these figures are stored in a specified
%directory
% if nargin==6
% plot_mc = 0;
% end
%%%
%----------------------read results  ------------------------%
fid = fopen(resultsfile,'r');
if fid == -1
    error(['Can not open result file: ' resultsfile]);
end

%get estimated global parameters from the first row of the result file
o = [];

while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'StartPos')),break,end
    % o
    result1 = regexp(tline,'LCR baseline shift:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        o = str2double(result1{1});
    end
end
%report errors if these values are not parsed successfully
if isempty(o)
    error(['Can not read estimated LCR baseline shift from ',resultsfile]);
end

%then read the results
results = textscan(fid,'%f %f %f %f %f %f %*f','treatAsEmpty', {'NA', 'na'});
fclose(fid);
chr_seg = results{1};
pstart_seg = results{2};
pend_seg = results{3};
cn_seg = results{4};
purity_seg = results{5};
score_seg = results{6};
clear results;

eval(['load ' Datafile]);

h=figure(1);
set(h,'visible','off');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 13 8])
FontSize = 17;
%--------------- plot figures ---------------------%
lcr_colors = [0.5 0.5 0.5;
             0 0.9 0;
             0 0 0.9;
             0.9 0 0];

for i = reshape(intersect(unique(data_chr_all),1:22),1,[]) %onlyl include autosome reshape(unique(data_chr_all),1,[])
    %process SNParray data
    tv = ismember(data_chr_all,i);
    data_lcr = data_lcr_all(tv);
    data_pos = data_spos_all(tv);
    min_pos = min(data_pos)-100;
    max_pos = max(data_pos)+100;
    
    indx1 = find(chr_seg==i);
    
    %plot
    clf;
    marker_size = 3;
    %---plot CN---
    subplot(3,1,1)
    hold on
    set(gca,'YGrid','on')
	set(gca,'FontSize',FontSize);
    %     axis ([-Inf Inf -0.05 7.05])
%     axis ([-Inf Inf -0.1 7.1])
    axis ([min_pos max_pos -0.1 7.1])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:1:7],'Box','on')
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('Copy number');
    for j=reshape(indx1,1,[])
        CN = cn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN > 7
            CN = 7;
        end
        line_style = 'r-';
        plot([data_pos(indx(1)) data_pos(indx(end))],[CN CN], line_style, 'LineWidth',3.0);
    end
    %replace '_' with '-' in barcode to display it correctly
	barcode_m = barcode;
    tmp = strfind(barcode_m,'_');
    barcode_m(tmp) = '-';
    set(gca,'XTick',[])
    title (['Chromosome ' num2str(i) ', ' barcode_m])
    
    subplot(3,1,2);
	set(gca,'FontSize',FontSize);
    hold on
    for j=reshape(indx1,1,[])
        CN = cn_seg(j);
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 1
            k = 1;
        else
            k = CN+1;
        end
        if k > 4
            k = 4;
        end
        plot(data_pos(tv),data_lcr(tv),'.','MarkerSize',marker_size, 'Color', lcr_colors(k,:));
    end
    % plot expected LCR mean values
    for j=reshape(indx1,1,[])
        CN = cn_seg(j);
        w = 1-purity_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN == 0
            CN = 0.001;
        end
        Y = w*2+(1-w)*CN;
        lcr_mean = log2(Y/2)+o;
        plot([data_pos(indx(1)) data_pos(indx(end))],[lcr_mean lcr_mean],'k-','LineWidth',1.5);
    end
    
    ylabel('LCR');
%     set(gca,'XTick',[]);
    set(gca,'Box','on')
    axis([min_pos max_pos -3 3])
    set(gca,'XTick',[])
%     axis([-Inf Inf -Inf Inf])

    %---plot Tumor fraction---
    subplot(3,1,3);
	set(gca,'FontSize',FontSize);
    hold on
    for j=reshape(indx1,1,[])
        line_style = 'r-';
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([data_pos(indx(1)) data_pos(indx(end))],[purity_seg(j) purity_seg(j)], line_style, 'LineWidth',1.5);
    end
    axis ([min_pos max_pos -0.03 1.03])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:0.5:1],'Box','on')
    ylabel('Cellularity');

%     subplot(4,1,4);
%     hold on
%     for j=reshape(indx1,1,[])
%         line_style = 'r-';
%         indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
%         if isempty(indx)
%             continue;
%         end
%         plot([data_pos(indx(1)) data_pos(indx(end))],[score_seg(j) score_seg(j)], line_style, 'LineWidth',1.5);
%     end
%     axis ([min_pos max_pos -10 110])
% %     set(gca,'XTick',[]);
%     set(gca,'YTick',[0:50:100],'Box','on')
%     ylabel('Score');

    %save figure
%     figpath = [plotsdir '\Chr_' num2str(i) '_' barcode];
    figpath = [plotsdir '/Chr_' num2str(i) '_' barcode '.png'];
    %     eval(['print -djpeg -r600 ' figpath ])
    eval(['print -dpng -r400 ' figpath ])

end