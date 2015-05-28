function [LL_all, num_snp] = CloneCNA_evaluate_BAF_LL(segments,beta)

global data_pos_sep
global data_bd_sep
global data_td_sep

global data_spos_ds_sep
global data_epos_ds_sep

LL_all = 0;
num_snp = 0;

cn = [2 0 1 3 4 5 6 7];

for i = 1:size(segments,1)
    chr_indx = segments(i,1);
    s_indx = segments(i,2);
    e_indx = segments(i,3);
    state_indx = segments(i,4);
    cp_indx = segments(i,5);
    St_pos = data_spos_ds_sep{chr_indx}(s_indx);
    Ed_pos = data_epos_ds_sep{chr_indx}(e_indx);
    
    snp_pos = data_pos_sep{chr_indx};
    tv = snp_pos >= St_pos & snp_pos <= Ed_pos;
    if sum(tv) == 0
        continue;
    end
    data_bd = data_bd_sep{chr_indx}(tv);
    data_td = data_td_sep{chr_indx}(tv);
    data_baf = data_bd./data_td;
    
    CN = cn(state_indx);
    if cp_indx == 0
        cp = 1;
    else
        cp = cp_indx;
    end
    majorCNs = ceil(CN/2):CN;
    if CN == 0
        Y = 0.001*beta(cp)+2*(1-beta(cp));
    else
        Y = CN*beta(cp)+2*(1-beta(cp));
    end
    Z = majorCNs'*beta(cp)+(1-beta(cp));
    centers = [Z/Y;1.0];% consider some purely homozygous SNPs
    het_wg = ones(length(centers),1);
    if CN > 0
        tv = centers == 0.5;
        het_wg(tv) = 1.3;
    end
    pie = [0.9*ones(length(centers)-1,1)/(length(centers)-1);0.1];
%     indi_update = [0 1];
%     LL = BMM(centers, pie, indi_update, het_wg, [data_bd data_td]);
    
    N = length(data_baf);
    mu_avg = mean(data_baf);
    sigma = sqrt((data_baf-mu_avg)'*(data_baf-mu_avg)/N)*ones(length(centers),1);
    indi_update = [0 1 1];
    if N == 1
        sigma(:) = 0.3;
        indi_update(3) = 0;
    end
    [LL, mu, pie, sigma] = GMM(centers, pie, sigma, indi_update, het_wg, data_baf);
%     disp(num2str(LL));
    
    if ~isnan(LL)
        num_snp = num_snp+length(data_baf);
        LL_all = LL_all+LL;
    end
    
%     figure(1);
%     cla();
%     minvalue = min(data_baf);
%     maxvalue = max(data_baf);
%     x1 = minvalue-1:0.1:maxvalue+1;
%     n = length(x1);
%     y1 = normpdf(repmat(x1',1,length(centers)),repmat(mu',n,1),repmat(sigma',n,1))*pie*0.5;
%     plot(x1,y1);
%     hold on;
    
%     x2 = minvalue:0.5:maxvalue;
%     [count,x2] = hist(data_baf,x2);
%     y2 = count/N;
%     bar(x2,y2);
    
    
end

end