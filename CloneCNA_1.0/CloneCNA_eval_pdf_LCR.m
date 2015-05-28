function results = CloneCNA_eval_pdf_LCR(data,mu,nu,sigma)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

temp = gammaln((nu+1)*0.5)-(nu+1)*0.5*log(1+(data-mu).^2/(nu*sigma^2))-gammaln(nu*0.5)-0.5*log(nu)-log(sigma)-0.5*log(pi);
results = exp(temp);
results(results<=eps) = eps;

end