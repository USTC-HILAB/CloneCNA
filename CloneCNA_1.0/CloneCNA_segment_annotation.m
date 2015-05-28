function [f_CN, f_cp] = CloneCNA_segment_annotation(data_lcr, data_baf, beta, nu, sigma, o)
% 09/15/2014 by Zhenhua
% This function is used to estimate copy number and clonal population for dismatched regions

% wg_prior = [1.3 1.7 1.7 1.5 1.5 1.1 1.1 1.1];
wg_prior = [1 1.7 1.7 1.4 1.3 1 0.8 0.7];
lcr_mean = median(data_lcr);
max_LL = [];
f_CN = -1;
f_cp = -1;
if length(beta) == 1
    f_CN = round((2^(lcr_mean-o+1)-2*(1-beta))/beta);
    f_cp = 1;
    return;
end
for cp = 1:length(beta)
    CN = round((2^(lcr_mean-o+1)-2*(1-beta(cp)))/beta(cp));
    if CN < 0
        CN = 0;
    end
    if CN == 0
        Y = 0.001*beta(cp)+2*(1-beta(cp));
    else
        Y = CN*beta(cp)+2*(1-beta(cp));  
    end
    majorCNs = ceil(CN/2):CN;
    Z = majorCNs'*beta(cp)+(1-beta(cp));
    
    mu_baf = Z/Y;
    
    if ~isempty(data_baf)        
        N = length(data_baf);
        mu_avg = mean(data_baf);
        centers = [mu_baf;1.0];
        pie = [0.9*ones(length(centers)-1,1)/(length(centers)-1);0.1];
        sigma_b = sqrt((data_baf-mu_avg)'*(data_baf-mu_avg)/N)*ones(length(centers),1);
        indi_update = [0 1 1];
        if N == 1
            sigma_b(:) = 0.3;
            indi_update(3) = 0;
        end
        het_wg = ones(length(centers),1);
        if CN > 0
            tv = centers == 0.5;
            het_wg(tv) = 1.5;
        end
        LL_baf = GMM(centers, pie, sigma_b, indi_update, het_wg, data_baf);
    else
        LL_baf = 0;
    end
    
    obslik_lcr = CloneCNA_eval_pdf_LCR(data_lcr,log2(Y/2)+o,nu,sigma);
    if CN > 7
        k = 8;
    else
        k = CN+1;
    end
    obslik_lcr = wg_prior(k)*obslik_lcr;
    LL_lcr = sum(log(obslik_lcr));
    LL = LL_baf+LL_lcr;
    if isempty(max_LL)
        max_LL = LL;
        f_CN = CN;
        f_cp = cp;
    elseif max_LL < LL
        max_LL = LL;
        f_CN = CN;
        f_cp = cp;
    end
end

end