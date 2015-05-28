function [p_states,num_loci,aCN,segments,alpha] = CloneCNA_process_results_new(beta,depend_table)
% 15/09/2014 by Zhenhua
%-----------------------------------------------------
%------overall information of the cancer sample------
%p_states: proportions of all hidden states
%num_loci: total number of loci investigated
%aCN: averaged copy number of clonal population
%segments: copy number segmentation results
%alpha: occurrence frequencies of clonal populations

global data_spos_ds_sep
global gamma_sep

tv_S = depend_table(:,2)~=0;
CN_mapping = depend_table(tv_S,3)'; % copy number of different entries

K = length(beta);

%initialize output parameters
num_loci = 0;
alpha = zeros(K+1,1);
segments = [];

%initialize intermediate variables
exp_num_states = [];
pos_dist = [];

for i = 1:length(gamma_sep) %for the ith chromosome
    post_probs = gamma_sep{i};
    data_pos = data_spos_ds_sep{i};
    
    %---handle p_states and num_loci---
    if isempty(exp_num_states) %initialization
        exp_num_states = zeros(size(post_probs,1),1);
    end
    exp_num_states = exp_num_states+sum(post_probs,2);
    num_loci = num_loci+size(post_probs,2);

    %---handle MAP states---
    %output predicted MAP states and subclonal populations 
    [temp,I] = max(post_probs,[],1);
    MAP_state = floor((I-1)/K)+1;
    MAP_sp = rem(I-1,K)+1;
    tv = MAP_state == max(MAP_state);
    MAP_state = MAP_state+1;
    MAP_state(tv) = 1;
    MAP_sp(tv) = 0;

    results = CloneCNA_segment_results(MAP_state,MAP_sp);
    segments = [segments; ones(size(results,1),1)*i results];    
    pos_dist = [pos_dist data_pos(results(:,2))-data_pos(results(:,1))+1]; 

    clear results;

end
pos_dist = pos_dist';

%---handle p_states---
p_states = zeros(length(CN_mapping),1);
for i = 1:length(CN_mapping)
    tv = segments(:,4) == i;
    if sum(tv) > 0
        p_states(i) = sum(pos_dist(tv))/sum(pos_dist);
    end
end

%---handle aCN---
aCN = zeros(1,K);
for k = 1:K
    pos_len = 0;
    w_aCN = 0;
    for i = 1:length(CN_mapping)
        tv = segments(:,4) == i & segments(:,5) == k;
        if sum(tv) > 0
            w_aCN = w_aCN+CN_mapping(i)*sum(pos_dist(tv));
            pos_len = pos_len+sum(pos_dist(tv));
        end
    end   
    aCN(k) = w_aCN/pos_len;
end
aCN = [aCN CN_mapping*p_states];

%---handle alpha---  
temp = sum(segments(:,3)-segments(:,2)+1);
for i = 1:K+1
    tv = segments(:,end) == i-1;
    alpha(i) = sum(segments(tv,3)-segments(tv,2)+1)/temp;
end

end