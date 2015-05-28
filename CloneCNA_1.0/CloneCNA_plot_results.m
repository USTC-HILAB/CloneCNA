function CloneCNA_plot_results(Datafile,resultsfile,plotsdir,barcode)
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

Chromosomes = intersect(unique(data_chr_all),1:22);
max_pos = zeros(1,length(Chromosomes));
max_lrc = max(data_lcr_all);
min_lcr = min(data_lcr_all);

for i = 1:length(Chromosomes)
    tv1 = data_chr_all == Chromosomes(i);
    max_pos(i) = max(data_spos_all(tv1));
end

ratio = max_pos/sum(max_pos);
xtick = cumsum([0 ratio(1:end-1)])+ratio/2;

clf;
line_style = 'k-';
LineWidth = 1.0;
MarkerSize = 4;
subplot(3,1,1);
hold on
set(gca,'YGrid','on');
set(gca,'FontSize',FontSize);
set(gca,'YTick',[0:1:7],'Box','on');
set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
ylabel('Copy number');
pre_x = 0;
chr_epos = zeros(length(Chromosomes),1);
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_spos_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN > 7
            CN = 7;
        end
        plot([x(indx(1)) x(indx(end))],[CN CN], 'r-', 'LineWidth',3.0);
    end
    chr_epos(i) = max(x);
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-0.1 7.1], line_style, 'LineWidth',LineWidth)
end
set(gca,'XTick',[])
axis([0 1 -0.1 7.1]);
barcode_m = barcode;
tmp = strfind(barcode_m,'_');
barcode_m(tmp) = '-';
title(barcode_m);

subplot(3,1,2);
set(gca,'FontSize',FontSize);
hold on
pre_x = 0;
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_spos_all(tv);
    data_lcr = data_lcr_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
    for j = reshape(indx1,1,[])
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
        plot(x(tv),data_lcr(tv),'.','MarkerSize',MarkerSize, 'Color', lcr_colors(k,:));
    end
    % plot expected LCR mean values
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        w = 1-purity_seg(j);
        if CN == 0
            CN = 0.001;
        end
        Y = w*2+(1-w)*CN;
        lcr_mean = log2(Y/2)+o;
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([x(indx(1)) x(indx(end))],[lcr_mean lcr_mean],'k-','LineWidth',1.5);
    end
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-3 3], line_style, 'LineWidth',LineWidth)
end
ylabel('LCR');
set(gca,'Box','on')
axis([0 1 -3 3])
% axis([0 1 -Inf Inf])
set(gca,'XTick',[])

subplot(3,1,3);
set(gca,'FontSize',FontSize);
hold on
pre_x = 0;
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_spos_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
    for j = reshape(indx1,1,[])
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([x(indx(1)) x(indx(end))],[purity_seg(j) purity_seg(j)], 'r-', 'LineWidth',1.5);
    end
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-0.03 1.03], line_style, 'LineWidth',LineWidth)
end
axis ([0 1 -0.03 1.03])
% set(gca,'XTick',[]);
set(gca,'YTick',[0:0.5:1],'Box','on')
ylabel('Cellularity');
xlabel('Chromosome');

set(gca,'XTick',xtick);
set(gca,'XTickLabel',mat2cell(Chromosomes',1,length(Chromosomes)));

% subplot(4,1,4);
% hold on
% pre_x = 0;
% for i = 1:length(Chromosomes)
%     tv = data_chr_all == Chromosomes(i);
%     data_pos = data_spos_all(tv);
%     x = data_pos*ratio(i)/max_pos(i)+pre_x;
%     indx1 = find(chr_seg == i);
%     for j = reshape(indx1,1,[])
%         indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
%         if isempty(indx)
%             continue;
%         end
%         plot([x(indx(1)) x(indx(end))],[score_seg(j) score_seg(j)], 'r-', 'LineWidth',1.5);
%     end
%     plot([max(x) max(x)], [-10 110], line_style, 'LineWidth',LineWidth)
%     pre_x = pre_x+ratio(i);  
% end
% for i = 1:length(Chromosomes)-1
%     plot([chr_epos(i) chr_epos(i)], [-10 110], line_style, 'LineWidth',LineWidth)
% end
% axis ([0 1 -10 110])
% set(gca,'YTick',[0:50:100],'Box','on')
% set(gca,'XTick',xtick);
% set(gca,'XTickLabel',mat2cell(Chromosomes',1,length(Chromosomes)));
% ylabel('Score');
% xlabel('Chromosome');

%save figure
figpath = [plotsdir '/' barcode '.png'];
eval(['print -dpng -r400 ' figpath ])