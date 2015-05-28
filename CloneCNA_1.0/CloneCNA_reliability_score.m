function score = CloneCNA_reliability_score(data_lcr,beta,nu,sigma,o,cn,cp)
% 15/09/2014 by Zhenhua
% This function is used to evaluate the reliability of CloneCNA results

if cn == 0
    cn = 0.001;
end

if cp == 0
    Y = 2;
else
    Y = cn*beta(cp)+2*(1-beta(cp));
end
mu_l = log2(Y/2)+o;

score_lcr = CloneCNA_eval_pdf_LCR(data_lcr,mu_l,nu,sigma)/CloneCNA_eval_pdf_LCR(mu_l,mu_l,nu,sigma);
score = median(score_lcr);

end
