function eta=soft_thresholding_C(x,lambda)
% complex soft-thresholding function

eta=(abs(x)> lambda).*(abs(x)-lambda).*(x)./abs(x+eps);

end