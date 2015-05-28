function [LL_all,CloneCNA_paras,p_states,num_loci,aCN,segments] = CloneCNA_screening(stepsize1,init_CloneCNA_paras, depend_table,thres1,max_iter1,verbose)

global data_lcr_sep
global data_spos_sep
global data_epos_sep

global data_lcr_ds_sep
global data_spos_ds_sep
global data_epos_ds_sep

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_CloneCNA_paras{1};
transmat_all = init_CloneCNA_paras{2};
beta_all = init_CloneCNA_paras{3};
nu_all = init_CloneCNA_paras{4};
sigma_all = init_CloneCNA_paras{5};
o_all = init_CloneCNA_paras{6};
indivec_all = init_CloneCNA_paras{7};

numex = length(data_lcr_sep); % each row is a sample
data_lcr_ds_sep = cell(1,numex);
data_spos_ds_sep = cell(1,numex);
data_epos_ds_sep = cell(1,numex);

for ex = 1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_lcr_sep{ex});
        data_lcr_ds_sep{ex} = data_lcr_sep{ex}(indx_ds);
        data_spos_ds_sep{ex} = data_spos_sep{ex}(indx_ds);
        data_epos_ds_sep{ex} = data_epos_sep{ex}(indx_ds);
    else %no ds
        data_lcr_ds_sep{ex} = data_lcr_sep{ex};
        data_spos_ds_sep{ex} = data_spos_sep{ex};
        data_epos_ds_sep{ex} = data_epos_sep{ex};
    end    
end

LL_all = [];
CloneCNA_paras = cell(1,7);
t_all = 0; 
if nargout >2
    p_states = [];
    num_loci = zeros(1,length(nu_all));
    aCN = cell(length(nu_all),1);
    segments = cell(length(nu_all),1);
end

for i=1:length(beta_all)
    if verbose
        tic
    end
    
%     disp(['parameter ' num2str(i)]);
    
    %1x1 cell
    init_CloneCNA_paras(1) = prior_all(i);
    init_CloneCNA_paras(2) = transmat_all(i);
    init_CloneCNA_paras(3) = beta_all(i);
    init_CloneCNA_paras(4) = nu_all(i);
    init_CloneCNA_paras(5) = sigma_all(i);
    init_CloneCNA_paras(6) = o_all(i);
    init_CloneCNA_paras(7) = indivec_all(i);
    
    [LL, prior, transmat, beta, nu, sigma, o, iterations] = CloneCNA_EM_Newton(init_CloneCNA_paras,depend_table,thres1,max_iter1,verbose);
        
    LL_all = [LL_all LL(end)];
    CloneCNA_paras{1} = [CloneCNA_paras{1} {prior}];
    CloneCNA_paras{2} = [CloneCNA_paras{2} {transmat}];
    CloneCNA_paras{3} = [CloneCNA_paras{3} {beta}];
    CloneCNA_paras{4} = [CloneCNA_paras{4} {nu}];
    CloneCNA_paras{5} = [CloneCNA_paras{5} {sigma}];
    CloneCNA_paras{6} = [CloneCNA_paras{6} {o}];
    CloneCNA_paras{7} = [CloneCNA_paras{7} init_CloneCNA_paras(7)];
    
    if nargout >2
        [temp,num_loci(i),aCN{i},segments{i}] = CloneCNA_process_results_new(beta,depend_table);
        p_states = [p_states temp];
    end

    if verbose
        t = toc;
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, iterations:' num2str(iterations) ', cp_num:' num2str(length(beta))]);
        disp(['beta:' num2str(reshape(beta,1,[]))]);
        disp(['nu:' num2str(nu) ', sigma:' num2str(sigma) ', o:' num2str(o)]);
        disp(['LL:' num2str(LL(end),'%5.1f') ', time:' num2str(t,'%5.1f')]);
        disp('--------------- screening report -----------------')
        t_all = t_all+t;
    end
    
end