function [CloneCNA_paras,best_indx] = CloneCNA_main(init_CloneCNA_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose,fn_nosuffix)
%--------------------------------------------------------------------%
%------------------>       version 1.0       <---------------------
%--------------------------------------------------------------------%
%03/10/2014 by Zhenhua
%CloneCNA main function,basically everything is done here
%--------------------------- screening -------------------------
global clamp_thres
global mc_w
global pre_best_paras
global NoSolutionFlag
clamp_thres = 1-1e-5;
mc_w = 0.8;
max_cp_num = 10;

pre_LL = -Inf;
pre_BAF_LL = -Inf;
pre_best_indx = [];
pre_CloneCNA_paras = [];
pre_best_paras = cell(1,6);
pre_aCN = [];

% bic_fid = fopen('INFO.bic', 'a+');
% fprintf(bic_fid,'----------------------------%s----------------------------\n',fn_nosuffix);
%-------------------------------------------------------------------
% Iteratively increase the number of clonal populations to find 
% the optimal model using Bayesian information criterion (BIC)
%-------------------------------------------------------------------
for cp_num = 1:max_cp_num
    cp_flag = 0;
    %-------------------------------------------------------------------
    %               ---> d-sampling screening <---
    %-------------------------------------------------------------------
    %initialize parameters
    thres_del = 0.009;

    init_CloneCNA_paras = CloneCNA_Init_paras(init_CloneCNA_paras,depend_table,cp_num,0);
    [LL,CloneCNA_paras,p_states,num_loci,aCN,segments] = CloneCNA_screening...
        (stepsize_ds,init_CloneCNA_paras,depend_table,thres_EM,max_iter,verbose);
       
    beta_all = cell2mat(CloneCNA_paras{3});
    nu_all = cell2mat(CloneCNA_paras{4});
    sigma_all = cell2mat(CloneCNA_paras{5});
    o_all = cell2mat(CloneCNA_paras{6});
    
%     CN = [0 1 2 3 4 5 6 7];
    N_genotype = [1 2 3 3 4 4 5 5];
    tv = depend_table(:,2) ~= 0;
    US_indx = depend_table(tv,1);
    
    DF_bias = (1-1./N_genotype(US_indx))*p_states;
    
    tv_del = depend_table(depend_table(:,2)~=0,3)<1;
    
    p_total_del = sum(p_states(tv_del,:),1);

    aCN = cell2mat(aCN)';
    if cp_num == 1
        tv = (p_total_del<thres_del) & (aCN(end,:) < 4.5);
    else
        tv = (p_total_del<thres_del) & (aCN(end,:) < 4.5) & ((aCN(end,:)-pre_aCN)/pre_aCN < 0.1);
    end
    if cp_num == 1
        NoSolutionFlag = false;
    end
    if ~any(tv)
        tv = ~tv;
        if cp_num == 1
            warning('Can not find a feasible solution with pre-defined criteria!');
            NoSolutionFlag = true;
        else
            break;
        end
    end
    candi_indx = find(tv);
    
    ratio = zeros(1,length(LL));
    % to select the optimal solution by investigating allelic read
    % depth data
    LL_baf = zeros(length(candi_indx),1);
    LL_all = zeros(length(candi_indx),1);
    num_snp = zeros(length(candi_indx),1);
    for i = 1:length(candi_indx)
        j = candi_indx(i);
        [LL_baf(i),num_snp(i)] = CloneCNA_evaluate_BAF_LL(segments{j},beta_all(:,j));
        LL_all(i) = LL(j)+LL_baf(i);
    end

    if cp_num == 1
        [temp,I] = max(LL_all);
        best_indx = candi_indx(I);
        candi_best = best_indx;
        
        thres1_ratio = 0.00005;
        %correction for possible signal noise
        for i = 1:length(candi_indx)
            j = candi_indx(i);
            temp1 = LL_all(I)-LL_all(i);
            temp2 = (log(num_loci(j)+num_snp(i))/2).*(DF_bias(best_indx)-DF_bias(j));
            if abs(temp1) <= 1 && abs(temp2) <= 0.01
                ratio(j) = 0;
            else
                ratio(j) = temp2/temp1;
            end
        end

        ds_candi1 = ratio >= thres1_ratio;
        if any(ds_candi1)
            [temp,I] = max(ratio(ds_candi1));
            indx = find(ds_candi1);
            best_indx = indx(I);
        end
    else
        tv1 = tv;
        if cp_num == 2
            tv1(candi_indx) = LL_baf >  pre_BAF_LL;
        end
        if sum(tv1) > 0
            [temp,I] = max(LL(tv1));
            candi_indx1 = find(tv1);
            best_indx = candi_indx1(I);
            candi_best = best_indx;
        else
            [temp,I] = max(LL(tv));
            best_indx = candi_indx(I);
            candi_best = best_indx;
        end
    end
    
%     % for test, save intermediate results
%     %-----------------------------------------------------------------------------------
%     fid = fopen('Model.details','a+');
%     fprintf(fid,'---------------------------%s----------------------------\n', fn_nosuffix);
%     fprintf(fid,'cp number:%d \n',cp_num);
%     for i = 1:size(aCN,2)
%         tv = i == candi_indx;
%         if i == best_indx
%             fprintf(fid,'2\t%f\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(end,i),ratio(i),LL(i),LL_baf(tv));
%         elseif i == candi_best
%             fprintf(fid,'1\t%f\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(end,i),ratio(i),LL(i),LL_baf(tv));
%         elseif sum(tv) > 0
%             fprintf(fid,'0\t%f\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(end,i),ratio(i),LL(i),LL_baf(tv));
%         else
%             fprintf(fid,'-1\t%f\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(end,i),ratio(i),LL(i),0);
%         end
%     end
% %     fprintf(fid,'\n');
%     fprintf(fid,'----------------------------------------------------------\n');
%     fclose(fid);
%     %-----------------------------------------------------------------------------------
    
    % Use Bayesian information criterion (BIC) to select the optimal model
    % With an increase in the number of subclonal populations (K to K+1), the number of parameters
    % that need to be estimated is increased by (2*K+1)*(S-1)^2+2*(S-1)+1
    
    alpha = 0.04;
    diff_th = 0;
    if cp_num == 1
        BIC_diff = -Inf;
        k = 0;
    else
        S = sum(depend_table(:,2)~=0);
        k = (2*(cp_num-1)+1)*(S-1)^2+2*(S-1)+1;
        BIC_diff = pre_LL-LL(best_indx)+alpha*0.5*k*log(num_loci(best_indx));
    end
%     fprintf(bic_fid,'%f\t%f\t%f\t%f\n',pre_LL,LL(best_indx),pre_LL-LL(best_indx),k*0.5*log(num_loci(best_indx)));
    if BIC_diff >= diff_th
        break;
    end
    
    tv = best_indx == candi_indx;
    if cp_num == 2 && LL_baf(tv) < pre_BAF_LL
        break;
    end
    pre_BAF_LL = LL_baf(tv);
    
%     
    fprintf(1, '\n');
    disp('--------------- new optimal parameters -----------------');
    disp(['cp number:' num2str(cp_num)]);
    disp(['beta:' num2str(reshape(beta_all(:,best_indx),1,[]))]);
    disp(['nu:' num2str(nu_all(best_indx)) ', sigma:' num2str(sigma_all(best_indx)) ', o:' num2str(o_all(best_indx))]);
    disp('--------------- new optimal parameters -----------------');
    fprintf(1, '\n');
    
    cp_flag = 1;
    pre_aCN = aCN(end,best_indx);
    pre_LL = LL(best_indx);
    pre_best_indx = best_indx;
    pre_CloneCNA_paras = CloneCNA_paras;
    init_CloneCNA_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]}];
    for i = 1:6
        pre_best_paras{i} = CloneCNA_paras{i}{best_indx};
    end
end
% fprintf(bic_fid,'----------------------------------------------------------\n\n');
% fclose(bic_fid);

% Now, use the optimal parameters to call CNA/LOH
%-------------------------------------------------------------------
if cp_flag == 1
    init_CloneCNA_paras = CloneCNA_Init_paras(pre_CloneCNA_paras,depend_table,cp_num,pre_best_indx);
else
    init_CloneCNA_paras = CloneCNA_Init_paras(pre_CloneCNA_paras,depend_table,cp_num-1,pre_best_indx);
end
[temp,CloneCNA_paras] = CloneCNA_screening(1,init_CloneCNA_paras,depend_table,5*thres_EM,20,verbose);
best_indx = 1;
%-------------------------------------------------------------------


function CloneCNA_paras = CloneCNA_Init_paras(init_CloneCNA_paras,depend_table,cp_num,best_indx)
%this function is used to initialize/process parameters for CloneCNA training,
%best_indx is used to indicate which parameter configuration (ususally
%multiple generated in previous screening procedure) are selected. If
%best_indx equals 0, parameters will be initialized
global clamp_thres
global var_l
global pre_best_paras

CloneCNA_paras = cell(1,7);
%parameter initialization
if best_indx == 0
    %---w---
    if isempty(init_CloneCNA_paras{3})
%         w_0 = linspace(0.001,0.9,10);
        if cp_num == 1
            beta_0 = [0.99 0.99 0.99 0.8 0.8 0.8 0.6 0.4 0.2 0.15 0.1];
            o_0 = [-0.6 -0.3 0 -0.6 -0.3 0 0 0 0 0 0];
        else
            pre_beta = pre_best_paras{3};
            if cp_num == 2
                ext_value1 = [0.1 0.15:0.1:pre_beta-0.05];
                ext_value2 = pre_beta+0.05:0.1:0.99;
                beta_0 = [ext_value1 repmat(pre_beta,1,length(ext_value2));repmat(pre_beta,1,length(ext_value1)) ext_value2];
            else                               
                beta1 = 0.1:0.05:pre_beta(1)-0.05;
                beta_0 = [beta1;repmat(pre_beta,1,length(beta1))];
                for i = 1:length(pre_beta)-1
                    beta1 = pre_beta(i)+0.05:0.05:pre_beta(i+1)-0.05;
                    temp = [repmat(pre_beta(1:i),1,length(beta1));beta1;repmat(pre_beta(i+1:end),1,length(beta1))];
                    beta_0 = [beta_0 temp];
                end
            end
            o_0 = repmat(pre_best_paras{6},1,size(beta_0,2));
        end
        if isempty(beta_0)
            beta_0 = repmat(linspace(0.1,0.99,cp_num)',1,3);
            o_0 = [-0.6 -0.3 0];
        end
    else
        beta_0 = init_CloneCNA_paras{3};
        o_0 = zeros(1,size(beta_0,2));
    end
    
    N = size(beta_0,2);
    CloneCNA_paras{3} = mat2cell(beta_0,size(beta_0,1),ones(1,N));
    
    %---nu--- 
    if isempty(init_CloneCNA_paras{4})
        if cp_num == 1
            nu_0 = 3;
        else
            nu_0 = pre_best_paras{4};
        end
    else
        nu_0 = init_CloneCNA_paras{4};
    end
    CloneCNA_paras{4} = repmat({nu_0},1,N);
    
    %---sigma--- 
    if isempty(init_CloneCNA_paras{5})
        if cp_num == 1
            sigma_0 = sqrt((nu_0-2)*var_l/nu_0);
        else
            sigma_0 = pre_best_paras{5};
        end
    else
        sigma_0 = init_CloneCNA_paras{5};
    end
    CloneCNA_paras{5} = repmat({sigma_0},1,N);
    
    %---o---
    CloneCNA_paras{6} = mat2cell(o_0,size(o_0,1),ones(1,N));

    %---pi---
    if isempty(init_CloneCNA_paras{1})
        S = sum(depend_table(:,2) ~= 0);
        K = size(beta_0,1);
        prior_0 = 1/((S-1)*K+1)*ones((S-1)*K+1,1);
    else
        prior_0 = init_CloneCNA_paras{1};
    end
    CloneCNA_paras{1} = repmat({prior_0},1,N);

    %---A---
    if isempty(init_CloneCNA_paras{2})
        S = sum(depend_table(:,2) ~= 0);
        K = size(beta_0,1);
        transmat_0 = norm_trans(ones((S-1)*K+1,(S-1)*K+1),clamp_thres);
    else
        transmat_0 = init_CloneCNA_paras{2};
    end
    CloneCNA_paras{2} = repmat({transmat_0},1,N);
    
    %---indicator vector---
    if isempty(init_CloneCNA_paras{7}) %indicator vector: '1' for update '0' fixed
        adj_all = ones(1,6);
    else
        adj_all = init_CloneCNA_paras{7};
    end
    CloneCNA_paras{7} = repmat({adj_all},1,N);
else %parse the results from previous screening
    for i = 1:length(best_indx)
        %--pi--
        CloneCNA_paras{1} = [CloneCNA_paras{1} init_CloneCNA_paras{1}(best_indx(i))];
        %--A--
        CloneCNA_paras{2} = [CloneCNA_paras{2} init_CloneCNA_paras{2}(best_indx(i))];
         %--beta--
        CloneCNA_paras{3} = [CloneCNA_paras{3} init_CloneCNA_paras{3}(best_indx(i))];
        %--nu--
        CloneCNA_paras{4} = [CloneCNA_paras{4} init_CloneCNA_paras{4}(best_indx(i))];
        %--sigma--
        CloneCNA_paras{5} = [CloneCNA_paras{5} init_CloneCNA_paras{5}(best_indx(i))];
        %--o--
        CloneCNA_paras{6} = [CloneCNA_paras{6} init_CloneCNA_paras{6}(best_indx(i))];
        %--indicator vector--
        CloneCNA_paras{7} = [CloneCNA_paras{7} init_CloneCNA_paras{7}(best_indx(i))];
    end
end
