function [d1, d2] = CalculateSoftThresholdDerivativeComplex(xc, lambda)

% function [dxR, dxI] = CalculateSoftThresholdDerivativeComplex(xc, lambda)
% This function calculates all the derivatives of the complex soft thresholding
% This is used to predict the correction term for the complex AMP
% algorithm.

xr = real(xc);
xi = imag(xc);
absx3over2 = (xr.^2+xi.^2).^(3/2)+eps;
indicatorabsx = (xr.^2+xi.^2>lambda^2);


dxR(:,1) = indicatorabsx.*(1- lambda*xi.^2./absx3over2);
dxR(:,2) = lambda*indicatorabsx.*xr.*xi./absx3over2;

%
dxI(:,1) = lambda*indicatorabsx.*xr.*xi./absx3over2;
dxI(:,2) =indicatorabsx.*(1- lambda*(xr.^2)./absx3over2); 

d1=dxR(:,1);
d2=dxI(:,2);

end