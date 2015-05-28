function [data_chr_all, data_spos_all, data_epos_all, data_lcr_all, data_chr_snp, data_pos_snp, data_bd_snp, data_td_snp] =...
    CloneCNA_load_preprocessData(countFile, refFile, depthFile, gcFile)
% 10/12/2014 by Zhenhua

%load read count data
fid = fopen(countFile, 'r');
if fid == -1
    error(['Can not open readcount file ' countFile]);
end
results = textscan(fid, '%f', 'treatAsEmpty', {'NA', 'na'}');
data_tc_all = results{1};
clear results;
fclose(fid);

fprintf(1,'Total %d targets are loaded from file "%s".\n',length(data_tc_all),countFile);
if length(data_tc_all) < 100000
    fprintf(1,'Warning: the number of targets loaded is too small, check the format of data file and whether data are completely loaded!\n');
end

fid = fopen(refFile, 'r');
if fid == -1
    error(['Can not open reference readcount file ' refFile]);
end
results = textscan(fid, '%f', 'treatAsEmpty', {'NA', 'na'}');
data_nc_all = results{1};
clear results;
fclose(fid);

fprintf(1,'Total %d targets are loaded from file "%s".\n',length(data_nc_all),refFile);
if length(data_nc_all) < 100000
    fprintf(1,'Warning: the number of targets loaded is too small, check the format of data file and whether data are completely loaded!\n');
end

%load GC-content data
fid = fopen(gcFile, 'r');
if fid == -1
    error(['Can not open GC-content file ' gcFile]);
end
results = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'treatAsEmpty', {'NA', 'na'}');
data_chr_all = results{1};
data_spos_all = results{2};
data_epos_all = results{3};
data_gc_all = results{4};
clear results;
fclose(fid);

tv = ismember(data_chr_all,{'X','x'});
data_chr_all(tv) = {'23'};
tv = ismember(data_chr_all,{'Y','y'});
data_chr_all(tv) = {'24'};
data_chr_all = str2double(data_chr_all);
clear tv;

% %load mappability data
% fid = fopen(mapFile, 'r');
% if fid == -1
%     error(['Can not open GC file ' mapFile]);
% end
% results = textscan(fid, '%*s%*f%*f%f', 'HeaderLines', 1, 'treatAsEmpty', {'NA', 'na'}');
% data_map_all = results{1};
% clear results;
% fclose(fid);

Chromosomes = intersect(unique(data_chr_all),1:22); % only use autosome
% tv = ismember(data_chr_all,Chromosomes) & ~(data_tc_all == 0 & data_nc_all == 0);
tv = ismember(data_chr_all,Chromosomes);
data_chr_all = data_chr_all(tv);
data_spos_all = data_spos_all(tv);
data_epos_all = data_epos_all(tv);
data_tc_all = data_tc_all(tv);
data_nc_all = data_nc_all(tv);
data_gc_all = data_gc_all(tv);
clear tv;

data_chr_all_n = [];
data_spos_all_n = [];
data_epos_all_n = [];
data_tc_all_n = [];
data_nc_all_n = [];
data_gc_all_n = [];

for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    tv1 = data_tc_all(tv) == 0 & data_nc_all(tv) == 0;
    if sum(tv1) == sum(tv)
        continue;
    end
    data_chr_all_n = [data_chr_all_n; data_chr_all(tv)];
    data_spos_all_n = [data_spos_all_n; data_spos_all(tv)];
    data_epos_all_n = [data_epos_all_n; data_epos_all(tv)];
    data_tc_all_n = [data_tc_all_n; data_tc_all(tv)];
    data_nc_all_n = [data_nc_all_n; data_nc_all(tv)];
    data_gc_all_n = [data_gc_all_n; data_gc_all(tv)];
end
data_chr_all = data_chr_all_n;
data_spos_all = data_spos_all_n;
data_epos_all = data_epos_all_n;
data_tc_all = data_tc_all_n;
data_nc_all = data_nc_all_n;
data_gc_all = data_gc_all_n;
clear data_chr_all_n data_spos_all_n data_epos_all_n data_tc_all_n data_nc_all_n data_gc_all_n;

Chromosomes = unique(data_chr_all);

%normalize for library size
data_tc_all = data_tc_all/sum(data_tc_all);
data_nc_all = data_nc_all/sum(data_nc_all);

data_ratio_all = (data_tc_all+eps)./(data_nc_all+eps);
tv = data_ratio_all >= 1/2^3 & data_ratio_all <= 2^3;
data_chr_all = data_chr_all(tv);
data_spos_all = data_spos_all(tv);
data_epos_all = data_epos_all(tv);
data_ratio_all = data_ratio_all(tv);
data_gc_all = data_gc_all(tv);

data_ratio_all = CloneCNA_ratio_correction(data_ratio_all, data_gc_all);

data_lcr_all = log2(data_ratio_all); % log count ratio

%load allelic read depth data
fid = fopen(depthFile, 'r');
if fid == -1
    error(['Can not open depth file ' depthFile]);
end
results = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'treatAsEmpty', {'NA', 'na'}');
data_chr_snp = results{1};
data_pos_snp = results{2};
data_bd_snp = results{3};
data_td_snp = results{4};
clear results;
fclose(fid);

tv = ismember(data_chr_snp,{'X','x'});
data_chr_snp(tv) = {'23'};
tv = ismember(data_chr_snp,{'Y','y'});
data_chr_snp(tv) = {'24'};
data_chr_snp = str2double(data_chr_snp);
clear tv;

fprintf(1,'Total %d SNPs are loaded from file "%s".\n',length(data_chr_snp),depthFile);
if length(data_chr_snp) < 10000
    fprintf(1,'Warning: the number of SNPs loaded is too small, check the format of data file and whether data are completely loaded!\n');
end

% filtering by read depth
minDepth = 10;
tv = ismember(data_chr_snp,Chromosomes) & data_td_snp >= minDepth;
data_chr_snp = data_chr_snp(tv);
data_pos_snp = data_pos_snp(tv);
data_bd_snp = data_bd_snp(tv);
data_td_snp = data_td_snp(tv);

% filter purely homozygous SNPs
data_baf_snp = data_bd_snp./data_td_snp;
tv = data_baf_snp > 0 & data_baf_snp < 1;
data_chr_snp = data_chr_snp(tv);
data_pos_snp = data_pos_snp(tv);
data_bd_snp = data_bd_snp(tv);
data_td_snp = data_td_snp(tv);

end

