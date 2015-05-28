function CloneCNA(TumorCountFile, NormalCountFile, TumorDepthFile, gcFile, outputDir, plotDir)
%--------------------------------------------------------------------%
%------------------>       version 1.0       <---------------------
%--------------------------------------------------------------------%
% 10/12/2014 by Zhenhua
% This is the first version of the CloneCNA method

%------Input and output------%
% TumorCountFile: tumor read count file
% NormalCountFile: normal or reference read count file
% TumorDepthFile: tumor allelic read depth file
% gcFile: GC-content file
% outputDir: directory of output files
% plotDir: directory of figure files

global current_version
global NoSolutionFlag
current_version = '1.0';

if nargin < 6
    error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m'] );
end

%parameters used in CloneCNA
%===============================================
%---for cn up to 7---
depend_table = [...
    %w1
    3 1 2;...
    1 1 0.001;...
    2 1 1;...
    4 1 3;...
    5 1 4;...
    6 1 5;...
    7 1 6;...
    8 1 7
    ];
%===============================================

thres_EM = 1e-4;
max_iter = 50;
verbose = 0;
init_CloneCNA_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]}]; % initial parameters will be assigned in the main function
                                                                   % parameters:pie,transmat,beta,nu,sigma,o
                                                                   
%initialization of global variable
global data_lcr_sep
global data_spos_sep
global data_epos_sep
global data_pos_sep
global data_bd_sep
global data_td_sep
global var_l

% record time
tic

max_copy = max(depend_table(:,3));

results = regexp(TumorCountFile,'/','split');
fn = results{end};
results = regexp(fn, '^(.+)\.+.+','tokens', 'once');
if isempty(results)
    fn_nosuffix = fn;
else
    fn_nosuffix = results{1};
    if ~isempty(strfind(fn_nosuffix,'.'))
        fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
    end
end

o_fid = fopen([outputDir '/' fn_nosuffix '.results'],'w');
if o_fid == -1
    error(['Can not open result file for ' fn_nosuffix]);
end

disp(['CloneCNA (version ' current_version ') is loading...'])

%--------------load and preprocess data--------------------
[data_chr_all, data_spos_all, data_epos_all, data_lcr_all, data_chr_snp, data_pos_snp, data_bd_snp, data_td_snp] = ...
    CloneCNA_load_preprocessData(TumorCountFile, NormalCountFile, TumorDepthFile, gcFile);

Chromosomes = reshape(unique(data_chr_all),1,[]);
            
nfilename = [outputDir '/' fn_nosuffix '_normalized'];

eval(['save ' nfilename '.mat data_chr_all data_spos_all data_epos_all data_lcr_all data_chr_snp data_pos_snp data_bd_snp data_td_snp']);

if size(data_spos_all,1)>size(data_spos_all,2) %make sure it's 1byN
    data_spos_all = data_spos_all';
end
if size(data_epos_all,1)>size(data_epos_all,2) %make sure it's 1byN
    data_epos_all = data_epos_all';
end
if size(data_lcr_all,1)>size(data_lcr_all,2) %make sure it's 1byN
    data_lcr_all = data_lcr_all';
end

%use at least 30000 data points for screening
stepsize_ds = max(floor(length(data_chr_all)/30000),1);

var_l = var(data_lcr_all);

%-------divide into different chromosomes----------
chr_num = length(Chromosomes);
data_spos_sep = cell(chr_num,1);
data_epos_sep = cell(chr_num,1);
data_lcr_sep = cell(chr_num,1);

data_pos_sep = cell(chr_num,1);
data_bd_sep = cell(chr_num,1);
data_td_sep = cell(chr_num,1);

for i = 1:chr_num
    tv = ismember(data_chr_all,Chromosomes(i));
    data_spos_sep{i} = data_spos_all(tv);
    data_epos_sep{i} = data_epos_all(tv);
    data_lcr_sep{i} = data_lcr_all(tv);
    tv = ismember(data_chr_snp,Chromosomes(i));
    data_pos_sep{i} = data_pos_snp(tv);
    data_bd = data_bd_snp(tv);
    data_td = data_td_snp(tv);
    data_baf = data_bd./data_td;               
    tv = data_baf < 0.5;           
    data_bd(tv) = data_td(tv)-data_bd(tv);   
    data_bd_sep{i} = data_bd;
    data_td_sep{i} = data_td;
end
clear data_chr_all data_spos_all data_epos_all data_lcr_all data_chr_snp data_pos_snp data_bd_snp data_td_snp;

%------------------ call CloneCNA --------------------
CloneCNA_paras = CloneCNA_main(init_CloneCNA_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose,fn_nosuffix);

beta = CloneCNA_paras{3}{1};
nu = CloneCNA_paras{4}{1};
sigma = CloneCNA_paras{5}{1};
o = CloneCNA_paras{6}{1};
   
[p_states,num_loci,aCN,segments,alpha] = CloneCNA_process_results_new(beta,depend_table);

cn_segs_all = zeros(size(segments,1),1); % copy number
cp_segs_all = zeros(size(segments,1),1); % clonal clusters
scores = zeros(size(segments,1),1); % reliability socre

cn = [2 0 1 3 4 5 6 7];
bp_len = 0;
cn_w = 0;
for i = 1:size(segments,1)
    chr_indx = segments(i,1);
    s_indx = segments(i,2);
    e_indx = segments(i,3);
    state_indx = segments(i,4);
    cp_indx = segments(i,5);
    St_pos = data_spos_sep{chr_indx}(s_indx);
    Ed_pos = data_epos_sep{chr_indx}(e_indx);
    
    snp_pos = data_pos_sep{chr_indx};
    tv = snp_pos >= St_pos & snp_pos <= Ed_pos;
    data_bd = data_bd_sep{chr_indx}(tv);
    data_td = data_td_sep{chr_indx}(tv);
    data_baf = data_bd./data_td;
    
    data_lcr = data_lcr_sep{chr_indx}(s_indx:e_indx);
    
    lcr_mean = median(data_lcr);
    if cp_indx == 0
        cp = 1;
    else
        cp = cp_indx;
    end
    
    if ~isempty(data_baf)
        majorCNs = ceil(cn(state_indx)/2):cn(state_indx);
        if cn(state_indx) == 0
            Y = 0.001*beta(cp)+2*(1-beta(cp));
        else
            Y = cn(state_indx)*beta(cp)+2*(1-beta(cp));
        end
        Z = majorCNs'*beta(cp)+(1-beta(cp));
        baf_mean = Z/Y;
        tv = data_baf < 0.97;
        tv = abs(median(data_baf(tv))-baf_mean) < 0.005;
        if sum(tv) == 0
            cn_est = 1;
        else
            cn_est = 0;
        end
    else
        cn_est = 0;
    end
    
    CN = round((2^(lcr_mean-o+1)-2*(1-beta(cp)))/beta(cp));
    
     % assign copy number and clonal population to each segment
    if CN ~= cn(state_indx) || cn_est == 1
        [CN,CP] = CloneCNA_segment_annotation(data_lcr, data_baf, beta, nu, sigma, o);
        cn_segs_all(i) = CN;
        cp_segs_all(i) = CP;
    else
        cn_segs_all(i) = cn(state_indx);
        cp_segs_all(i) = cp_indx;
    end
    
%     cn_segs_all(i) = cn(state_indx);
%     cp_segs_all(i) = cp_indx;
    if cn_segs_all(i) == 2
        cp_segs_all(i) = 0;
    end

    bp_len = bp_len+(Ed_pos-St_pos+1);   
    cn_w = cn_w+(Ed_pos-St_pos+1)*cn_segs_all(i);
    scores(i) = CloneCNA_reliability_score(data_lcr,beta,nu,sigma,o,cn_segs_all(i),cp_segs_all(i));
end

seg_len = segments(:,3)-segments(:,2)+1;
m_score = sum(scores.*seg_len)/sum(seg_len);
max_score = max(scores);

scores = scores*100*2/(max_score+m_score);
scores(scores>100) = 100;

ACN = cn_w/bp_len;

listfile = 'LOG.txt';
fp_list = fopen(listfile,'a+');
if fp_list == -1
    warning('Can not open the file: "%s"!\n',listfile);
else
    fprintf(fp_list,'Version:\tDate\tTime\tSample\tTumor purity\tCellularity of subclonal populations\tTumor ploidy\tnu\tSigma\to\n');
    fprintf(fp_list,'CloneCNA%s\t%s\t%s\t%s\t%f\t',current_version,datestr(clock,'mmm-dd-yyyy'),datestr(clock,'HH:MM:SS'),fn_nosuffix,beta(end));
    for i = 1:length(beta)
        if i == length(beta)
            fprintf(fp_list,'%f\t',beta(i));
        else
            fprintf(fp_list,'%f,',beta(i));
        end
    end
    fprintf(fp_list,'%f\t%d\t%f\t%f\n',ACN,nu,sigma,o);
end
fclose(fp_list);

%-------------- output summary of the results--------------
fprintf(o_fid,'---------------------------------------------------------------\n');
fprintf(o_fid,['             Summary of CloneCNA results (version ' current_version ')          \n']);
if NoSolutionFlag
    fprintf(o_fid,'Warning: Prediction results may be inaccurate due to the failure\n');
    fprintf(o_fid,'in finding optimal initial global parameters!\n');
end

fprintf(o_fid,'General information of this cancer sample:                      \n');
fprintf(o_fid,'   Proportion of abnormal cells in the sample: %6.4f\n',beta(end));
fprintf(o_fid,'   Average copy number: %1.2f\n',ACN);
fprintf(o_fid,'   LCR baseline shift: %6.2f\n',o);
fprintf(o_fid,'   Sigma of t-distribution: %6.4f\n',sigma);
fprintf(o_fid,'   Degree of freedom of t-distribution: %d\n',nu);

fprintf(o_fid,'   Cellularity of subclonal populations in the sample:');
for i = 1:length(beta)
    fprintf(o_fid,' %6.4f',beta(i));
end
fprintf(o_fid,'\n');

fprintf(o_fid,'---------------------------------------------------------------\n');
fprintf(o_fid,'\n');

%--------------output state assignment in segments--------------
fprintf(o_fid,'Chr\tStartPos\tEndPos\tCN\tCellularity\tScore\n');
for i = 1:size(segments,1)
    chr_indx = segments(i,1);
    s_indx = segments(i,2);
    e_indx = segments(i,3);

    s_pos = data_spos_sep{chr_indx}(s_indx);
    e_pos = data_epos_sep{chr_indx}(e_indx);
    if cp_segs_all(i) == 0
        cellularity = -1;
    else
        cellularity = beta(cp_segs_all(i));
    end
    fprintf(o_fid,'%d\t%d\t%d\t%d\t%f\t%3.1f\n',Chromosomes(chr_indx),...
                s_pos,e_pos,cn_segs_all(i),cellularity,scores(i));
end
fclose(o_fid);

disp ('-----Plot CloneCNA results now-----');

Datafile = [outputDir '/' fn_nosuffix '_normalized.mat'];
resultsfile = [outputDir '/' fn_nosuffix '.results'];  
plotsdir = [plotDir '/' fn_nosuffix];
s = mkdir(plotsdir);
if ~s %failed to make a directory
    error(['Can not make directory: ' plotsdir]);
else
    CloneCNA_plot_normalized_results(Datafile,resultsfile,plotsdir,fn_nosuffix);
    CloneCNA_plot_results(Datafile,resultsfile,plotsdir,fn_nosuffix);
end

t = toc;
m = floor(t/60);
s = t-m*60;
disp (['sample ' fn_nosuffix ' is done, time used: ' num2str(m) ' minutes and ' num2str(s) 'seconds']);
clear all
close all

end