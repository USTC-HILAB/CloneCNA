function [LL, mu, pie] = BMM(centers, pie, indi_update, data)

zero_final = eps;
th = 1.0e-5;

% initialize
N = size(data,1);
% pie = repmat(1/center_num,center_num,1);
mu = centers;
center_num = length(centers);

if size(mu,1) < size(mu,2)
    mu = mu';
end
if size(pie,1) < size(pie,2)
    pie = pie';
end

iter_num = 0;
max_iter = 100;
LL_all = [];
pre_LL = -inf;
% tic
% EM Algorithm
while 1
    % expectation step
    pdf_values = binopdf(repmat(data(:,1),1,center_num),repmat(data(:,2),1,center_num),repmat(mu',N,1));
    pdf_values(pdf_values < zero_final) = zero_final;
    numerator = repmat(pie',N,1).*pdf_values;
    denominator = repmat(sum(numerator,2),1,center_num);
    gamma = numerator./denominator;
    LL = sum(log(sum(numerator,2)));
    
    delta_LL = abs(LL-pre_LL);
    avg_LL = (abs(LL)+abs(pre_LL)+eps)/2;
    if delta_LL/avg_LL < th
        break;
    end
    
    pre_LL = LL;
    LL_all = [LL_all LL];
    % maximization step
    if indi_update(1) == 1
        mu = gamma'*(data(:,1)./data(:,2))./sum(gamma)';
    end
    if indi_update(2) == 1
        pie = (sum(gamma)/N)';
    end
    iter_num = iter_num+1;
    
%     disp('pie: ');
%     for i = 1:center_num
%         disp(num2str(pie(i)));
%     end
%     disp('mu: ');
%     for i = 1:center_num
%         disp(num2str(mu(i)));
%     end
    
    if iter_num >= max_iter
        break;
    end    
end
% toc

% disp(['log-likelihood: ' num2str(LL)]);
% disp('mu: ');
% for i = 1:center_num
%     disp(num2str(mu(i)));
% end
% disp('pie: ');
% for i = 1:center_num
%     disp(num2str(pie(i)));
% end

end
