function [LL, mu, pie, sigma] = GMM(centers, pie, sigma, indi_update, het_wg, data)

zero_final = eps;
th = 1.0e-5;

% initialize means and sigmas
N = length(data);
center_num = length(centers);
% mu_avg = mean(data);
% n = rem(N,center_num);
% temp = sort(data);
% mu = mean(reshape(temp(n+1:end),[],center_num))';
% sigma = repmat(sqrt((data-mu_avg)'*(data-mu_avg)/N),center_num,1);
% pie = repmat(1/center_num,center_num,1);

mu = centers;
if size(mu,1) < size(mu,2)
    mu = mu';
end
if size(pie,1) < size(pie,2)
    pie = pie';
end
if size(sigma,1) < size(sigma,2)
    sigma = sigma';
end
if size(het_wg,1) < size(het_wg,2)
    het_wg = het_wg';
end

iter_num = 0;
max_iter = 100;
LL_all = [];
pre_LL = -inf;
% tic
% EM Algorithm
while 1
    % expectation step
    pdf_values = normpdf(repmat(data,1,center_num),repmat(mu',N,1),repmat(sigma',N,1)).*repmat(het_wg',N,1);
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
        mu = gamma'*data./sum(gamma)';
        mu(end) = 1.0;
    end
    if indi_update(3) == 1
        temp = repmat(data,1,center_num) - repmat(mu',N,1);
        sigma = sqrt(diag(gamma'*(temp.^2))./sum(gamma)');
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
% disp('pie: ');
% for i = 1:center_num
%     disp(num2str(pie(i)));
% end
% disp('mu: ');
% for i = 1:center_num
%     disp(num2str(mu(i)));
% end
% disp('sigma: ');
% for i = 1:center_num
%     disp(num2str(sigma(i)));
% end

% figure(1);
% cla();
% subplot(2,1,1);
% minvalue = min(data);
% maxvalue = max(data);
% x1 = minvalue-1:0.1:maxvalue+1;
% n = length(x1);
% y1 = normpdf(repmat(x1',1,center_num),repmat(mu',n,1),repmat(sigma',n,1))*pie*0.5;
% plot(x1,y1);
% hold on;
% 
% x2 = minvalue:0.5:maxvalue;
% [count,x2] = hist(data,x2);
% y2 = count/N;
% bar(x2,y2);
% 
% 
% subplot(2,1,2);
% plot(LL_all);
% xlabel('iteration');
% ylabel('log likelihood');

end
