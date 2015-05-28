function [obslik,condi_probs_fluct] = CloneCNA_get_obslik(data_lcr,beta,nu,sigma,o,depend_table)

N = length(data_lcr); %number of data points

ns = 2; %copy number of stromal cells
tv_S = depend_table(:,2)~=0;
Nc = depend_table(tv_S,3); %vector of copy numbers of different entries

temp1 = repmat(Nc,1,length(beta));
temp2 = repmat(beta',length(Nc),1);
Y = temp1.*temp2+ns*(1-temp2);
mu_l = log2(Y/2)+o;

S = sum(tv_S);
K = length(beta);

obslik = zeros((S-1)*K+1,N);
condi_probs_fluct = zeros((S-1)*K+1,N);

% copy neutral state
obslik_lcr = CloneCNA_eval_pdf_LCR(data_lcr,o,nu,sigma);

fluct_prob = 0.001;
obslik((S-1)*K+1,:) = (1-fluct_prob)*obslik_lcr+fluct_prob/6;
condi_probs_fluct((S-1)*K+1,:) = (fluct_prob/6)./obslik((S-1)*K+1,:);

for j = 1:K
    for i = 2:size(Y,1)
        obslik_lcr = CloneCNA_eval_pdf_LCR(data_lcr,mu_l(i,j),nu,sigma);      
%         if Nc(i) < 1
%             obslik_lcr = 0.9*obslik_lcr;
%         end
        obslik((i-2)*K+j,:) = (1-fluct_prob)*obslik_lcr+fluct_prob/6;
        condi_probs_fluct((i-2)*K+j,:) = (fluct_prob/6)./obslik((i-2)*K+j,:);
    end
end

end

