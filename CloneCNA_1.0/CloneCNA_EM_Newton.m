function [LL, prior, transmat, beta, nu, sigma, o, nrIterations] = ...
    CloneCNA_EM_Newton(init_CloneCNA_paras,depend_table, thresh, max_iter,verbose)
% 15/09/2014 by Zhenhua
% This new EM algorithm uses univariate method to update 
% parameters sperately, Newton¨CRaphson method is used
global clamp_thres

previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];

prior = init_CloneCNA_paras{1};
transmat = init_CloneCNA_paras{2};
beta = init_CloneCNA_paras{3};
nu = init_CloneCNA_paras{4};
sigma = init_CloneCNA_paras{5};
o = init_CloneCNA_paras{6};

while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [loglik, exp_num_trans, exp_num_visits1,beta_u,nu_u,sigma_u,o_u] = ...
        CloneCNA_compute_ess(prior,transmat,beta,nu,sigma,o,depend_table,verbose);
    
    [converged, decrease] = em_converged_m(loglik, previous_loglik, verbose,thresh);
    if decrease
%         converged = 1;
%         continue;
    end
    
    % update parameters
    if init_CloneCNA_paras{7}(1)
        prior = norm_trans(exp_num_visits1',0)';
    end
    if init_CloneCNA_paras{7}(2) && ~isempty(exp_num_trans)
        % clamp_thres = 1-1e-4;
        transmat = norm_trans(exp_num_trans,clamp_thres);
    end
    if init_CloneCNA_paras{7}(3) %update w here
        beta = beta_u;
    end
    if init_CloneCNA_paras{7}(4) %update nu here
        nu = nu_u;
    end
    if init_CloneCNA_paras{7}(5) %update sigma here
        sigma = sigma_u;
    end
    if init_CloneCNA_paras{7}(6) %update sigma here
        o = o_u;
    end
    
    if verbose
        disp(['beta:' num2str(reshape(beta,1,[]))]);
        disp(['nu:' num2str(nu) ', sigma:' num2str(sigma) ', o:' num2str(o)]);
        fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik);
    end
    
    num_iter =  num_iter + 1;
    previous_loglik = loglik;
    LL = [LL loglik];
end
nrIterations = num_iter - 1;

end

%--------------------------------------------------------------------------
function [loglik, exp_num_trans, exp_num_visits1,beta_u,nu_u,sigma_u,o_u] = ...
    CloneCNA_compute_ess(prior,transmat,beta,nu,sigma,o,depend_table,verbose)
global data_lcr_ds_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lcr_ds_sep); % each row is a sample
S_all = length(transmat); % number of all states 
exp_num_trans = zeros(S_all,S_all);
exp_num_visits1 = zeros(S_all,1);

%-----------------------E step-----------------------------
gamma_sep = cell(1,numex);
condi_probs_fluct_sep = cell(1,numex);
loglik = 0;
N = 0; % the size of the whole data set

for ex=1:numex %
    
    % conditional probabilities
    [obslik,condi_probs_fluct] = CloneCNA_get_obslik(data_lcr_ds_sep{ex},beta,nu,sigma,o,depend_table);
    % Forward and Backward algorithm
    [temp1, gamma, current_ll, temp2, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    
    clear temp1 temp2;
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    
    gamma_sep{ex} = gamma;
    clear gamma;
    condi_probs_fluct_sep{ex} = condi_probs_fluct;
    clear condi_probs_fluct;
end

%-----------------------M step-----------------------------
%update global parameters
max_iter = 10;

%update beta
beta_tol = 0.01;
beta_u = CloneCNA_update_beta(beta,nu,sigma,o,depend_table,beta_tol,max_iter);
% beta_u = beta;

%update sigma
sigma_tol = 0.01;
sigma_u = CloneCNA_update_sigma(beta_u,nu,sigma,o,depend_table,sigma_tol,max_iter);

%update o
o_tol = 0.01;
o_u = CloneCNA_update_o(beta_u,nu,sigma_u,o,depend_table,o_tol,max_iter);
% o_u = 0;

%update nu
nu_tol = 0.5;
nu_u = CloneCNA_update_nu(beta_u,nu,sigma_u,o_u,depend_table,nu_tol,max_iter);

end

%--------------------------------------------------------------------------
function beta_u = CloneCNA_update_beta(beta,nu,sigma,o,depend_table,beta_tol,max_iter)
global data_lcr_ds_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lcr_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 8; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %row vector of copy numbers of different entries

normal_indx = find(Nc == 2);

iter = 0;
while 1
    %first order differential
    ELL_L_D_1 = zeros(length(beta),1);
    %second order differential
    ELL_L_D_2 = zeros(length(beta),1);
    
    temp1 = repmat(Nc,1,K);
    temp2 = repmat(beta',length(Nc),1);
    Y = temp1.*temp2+ns*(1-temp2);
    
    mu_l = log2(Y/2)+o;
    mu_l_D1 = repmat(Nc-ns,1,K)./(Y*log(2));
    mu_l_D2 = repmat(-(Nc-ns).^2,1,K)./(Y.^2*log(2));

    for ex=1:numex %
        obs_lcr = data_lcr_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)              
            if i == normal_indx
                continue;
            end
            for j = 1:size(Y,2)
                k = (i-2)*K+j;
                ELL_L_D_1(j) = ELL_L_D_1(j)+post_probs(k,:)*((nu+1)*(obs_lcr-mu_l(i,j))*mu_l_D1(i,j)./(sigma^2*nu+(obs_lcr-mu_l(i,j)).^2))';               
                ELL_L_D_2(j) = ELL_L_D_2(j)+post_probs(k,:)*(nu+1)*((2*(obs_lcr-mu_l(i,j)).^2*mu_l_D1(i,j)^2./sigma^2-...
                    (mu_l_D1(i,j)^2+(obs_lcr-mu_l(i,j))*mu_l_D2(i,j)).*(nu+(obs_lcr-mu_l(i,j)).^2/sigma^2))./(sigma*nu+(obs_lcr-mu_l(i,j)).^2/sigma).^2)';                   
            end
        end
    end
    beta_adj = -ELL_L_D_1./ELL_L_D_2;
    
    % now determine if the update violates the constrains
    % let beta_u = beta+c*beta_adj, c is the minimal coefficient from a list of
    % coefficients c_all that activiate but do not violate the constrains
    % 0.01<=beta(1)<beta(2)<...<beta(K)<=1
    c_p = 1;
    
    %-beta(1)<=-0.01 => beta(1)+c*beta_adj(1)>=0.01
    if -beta_adj(1) > 0 % not a feasible direction
        temp = (0.01-beta(1))/beta_adj(1);
        if temp < 1
%             disp(['Constrain: beta(' num2str(1) ')>' num2str(0) ') is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    for j = 2:K
        %-beta(j)+beta(j-1)<=-0.05 => beta(j)+c*beta_adj(j)>=beta(j-1)+c*beta_adj(j-1)+0.05
        if -beta_adj(j)+beta_adj(j-1) > 0 % not a feasible direction
            temp = (beta(j)-beta(j-1)-0.05)/(beta_adj(j-1)-beta_adj(j));
            if temp < 1
%                 disp(['Constrain: beta(' num2str(j-1) ')<beta(' num2str(j) ') is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
    end
    %beta(K)-0.99<=0 => beta(K)+c*beta_adj(K)<=0.99
    if beta_adj(K) > 0 % not a feasible direction
        temp = (0.99-beta(K))/beta_adj(K);
        if temp < 1
%             disp(['Constrain: beta(' num2str(K) ')<=1 is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    beta_u = beta+c_p*beta_adj;
    tv = isnan(beta_u);
    beta_u(tv) = beta(tv);
    
    iter = iter+1;
    
%     if sum(unConverged) == 0 || iter > max_iter
    if max(abs(beta_u-beta)) < beta_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        beta = beta_u;
    end
end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function nu_u = CloneCNA_update_nu(beta,nu,sigma,o,depend_table,nu_tol,max_iter)

global data_lcr_ds_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lcr_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 8; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries

normal_indx = find(Nc == 2);

temp1 = repmat(Nc,1,K);
temp2 = repmat(beta',length(Nc),1);
Y = temp1.*temp2+ns*(1-temp2);
clear temp1 temp2;

mu_l = log2(Y/2)+o;

iter = 0;
while 1
    %first order differential
    ELL_L_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    
    for ex = 1:numex
        obs_lcr = data_lcr_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)              
            if i == normal_indx
                ELL_L_D_1 = ELL_L_D_1+post_probs(end,:)*0.5*(psi((nu+1)*0.5)-log(1+(obs_lcr-o).^2/(nu*sigma^2))+(nu+1)*(obs_lcr-o).^2./((nu*sigma)^2+...
                         nu*(obs_lcr-o).^2)-psi(nu*0.5)-1/nu)';            
                ELL_L_D_2 = ELL_L_D_2+post_probs(end,:)*0.5*(0.5*psi(1,(nu+1)*0.5)+(obs_lcr-o).^2./(nu^2*sigma^2+nu*(obs_lcr-o).^2)+...
                         (obs_lcr-o).^2.*(nu^2+nu*(obs_lcr-o).^2/sigma^2-(nu+1)*(2*nu+(obs_lcr-o).^2/sigma^2))./(nu^2*sigma+nu*(obs_lcr-o).^2/sigma).^2-...
                         0.5*psi(1,nu*0.5)+1/nu^2)';  
                continue;
            end
            for j = 1:size(Y,2)
                k = (i-2)*K+j;
                ELL_L_D_1 = ELL_L_D_1+post_probs(k,:)*0.5*(psi((nu+1)*0.5)-log(1+(obs_lcr-mu_l(i,j)).^2/(nu*sigma^2))+(nu+1)*(obs_lcr-mu_l(i,j)).^2./((nu*sigma)^2+...
                         nu*(obs_lcr-mu_l(i,j)).^2)-psi(nu*0.5)-1/nu)';            
                ELL_L_D_2 = ELL_L_D_2+post_probs(k,:)*0.5*(0.5*psi(1,(nu+1)*0.5)+(obs_lcr-mu_l(i,j)).^2./(nu^2*sigma^2+nu*(obs_lcr-mu_l(i,j)).^2)+...
                         (obs_lcr-mu_l(i,j)).^2.*(nu^2+nu*(obs_lcr-mu_l(i,j)).^2/sigma^2-(nu+1)*(2*nu+(obs_lcr-mu_l(i,j)).^2/sigma^2))./(nu^2*sigma+nu*(obs_lcr-mu_l(i,j)).^2/sigma).^2-...
                         0.5*psi(1,nu*0.5)+1/nu^2)';          
            end           
        end
    end
    
    nu_adj = -ELL_L_D_1/ELL_L_D_2;
    
    c_p = 1;
    
    %-nu<=-3 => nu+c*nu_adj>=3
    if -nu_adj > 0 % not a feasible direction
        temp = (3-nu)/nu_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    %nu<=10000 => nu+c*nu_adj<=10000
    if nu_adj > 0 % not a feasible direction
        temp = (10000-nu)/nu_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    nu_u = nu+c_p*nu_adj;
    if isnan(nu_u)
        nu_u = nu;
    end
    
    iter = iter+1;
    
    if abs(nu_u-nu) < nu_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        nu_u = round(nu_u);
        break;
    else
        nu = round(nu_u);
    end
end

end

%--------------------------------------------------------------------------
function sigma_u = CloneCNA_update_sigma(beta,nu,sigma,o,depend_table,sigma_tol,max_iter)

global data_lcr_ds_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lcr_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 8; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries

normal_indx = find(Nc == 2);

temp1 = repmat(Nc,1,K);
temp2 = repmat(beta',length(Nc),1);
Y = temp1.*temp2+ns*(1-temp2);
clear temp1 temp2;

mu_l = log2(Y/2)+o;

iter = 0;
while 1
    %first order differential
    ELL_L_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    
    for ex = 1:numex
        obs_lcr = data_lcr_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)              
            if i == normal_indx
                ELL_L_D_1 = ELL_L_D_1+post_probs(end,:)*(nu*(-sigma^2+(obs_lcr-o).^2)./(sigma*(nu*sigma^2+(obs_lcr-o).^2)))';         
                ELL_L_D_2 = ELL_L_D_2+post_probs(end,:)*nu*((sigma^2-(obs_lcr-o).^2)./(sigma^2*(nu*sigma^2+(obs_lcr-o).^2))-2*(obs_lcr-o).^2*(1+nu)./(nu*sigma^2+(obs_lcr-o).^2).^2)'; 
                continue;
            end
            for j = 1:size(Y,2)
                k = (i-2)*K+j;
                ELL_L_D_1 = ELL_L_D_1+post_probs(k,:)*(nu*(-sigma^2+(obs_lcr-mu_l(i,j)).^2)./(sigma*(nu*sigma^2+(obs_lcr-mu_l(i,j)).^2)))';         
                ELL_L_D_2 = ELL_L_D_2+post_probs(k,:)*nu*((sigma^2-(obs_lcr-mu_l(i,j)).^2)./(sigma^2*(nu*sigma^2+(obs_lcr-mu_l(i,j)).^2))-2*(obs_lcr-mu_l(i,j)).^2*(1+nu)./(nu*sigma^2+(obs_lcr-mu_l(i,j)).^2).^2)';          
            end           
        end
    end
    
    sigma_adj = -ELL_L_D_1/ELL_L_D_2;
    
    c_p = 1;
    
    %-sigma<=-0.01 => sigma+c*sigma_adj>=0.01
    if -sigma_adj > 0 % not a feasible direction
        temp = (0.01-sigma)/sigma_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    %sigma<=0.5 => sigma+c*sigma_adj<=0.5
    if sigma_adj > 0 % not a feasible direction
        temp = (0.5-sigma)/sigma_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    sigma_u = sigma+c_p*sigma_adj;
    if isnan(sigma_u)
        sigma_u = sigma;
    end
    
    iter = iter+1;
    
    if abs(sigma_u-sigma) < sigma_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        sigma = sigma_u;
    end
end

end

%--------------------------------------------------------------------------
function o_u = CloneCNA_update_o(beta,nu,sigma,o,depend_table,o_tol,max_iter)
global data_lcr_ds_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lcr_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 8; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries

normal_indx = find(Nc == 2);

temp1 = repmat(Nc,1,K);
temp2 = repmat(beta',length(Nc),1);
Y = temp1.*temp2+ns*(1-temp2);
clear temp1 temp2;

iter = 0;
while 1
    %first order differential
    ELL_L_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    
    mu_l = log2(Y/2)+o;
    
    for ex = 1:numex
        obs_lcr = data_lcr_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)              
            if i == normal_indx
                ELL_L_D_1 = ELL_L_D_1+post_probs(end,:)*((nu+1)*(obs_lcr-o)./(sigma^2*nu+(obs_lcr-o).^2))';  
                ELL_L_D_2 = ELL_L_D_2+post_probs(end,:)*(nu+1)*((2*(obs_lcr-o).^2/sigma^2-...
                    (nu+(obs_lcr-o).^2/sigma^2))./(sigma*nu+(obs_lcr-o).^2/sigma).^2)'; 
                continue;
            end
            for j = 1:size(Y,2)
                k = (i-2)*K+j;
                ELL_L_D_1 = ELL_L_D_1+post_probs(k,:)*((nu+1)*(obs_lcr-mu_l(i,j))./(sigma^2*nu+(obs_lcr-mu_l(i,j)).^2))';  
                ELL_L_D_2 = ELL_L_D_2+post_probs(k,:)*(nu+1)*((2*(obs_lcr-mu_l(i,j)).^2/sigma^2-...
                    (nu+(obs_lcr-mu_l(i,j)).^2/sigma^2))./(sigma*nu+(obs_lcr-mu_l(i,j)).^2/sigma).^2)'; 
            end           
        end
    end
    
    o_u = o-ELL_L_D_1/ELL_L_D_2;
    if isnan(o_u)
        o_u = o;
    end
    
    iter = iter+1;
    
    if abs(o_u-o) < o_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        o = o_u;
    end
end

end

