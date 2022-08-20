function out1 = genericAMP( y,A,colnormA,Eta,Etader, par )
% This algorithm is an extension of AMP and SE to the complex setting
% Inputs:
%   y :          observations
%   A :          a function handle that represents matrix A
%   Eta :        a function handle which is a generic denoiser, 
%                xhat=Eta(temp,sigma) 
%   Etader :     a function handle which is the derivative function of the 
%                denoise function Eta. This function takes a vector and a 
%                value of thresholding (similar as the parameters that the 
%                function Eta needs), and should return two vectors, the 
%                first vector is the partial derivative of real(Eta(temp,sigma))
%                with respect to real(temp), the second vector is the partial
%                derivative of imag(Eta(temp,sigma)) with respect to imag(temp).
%                If you can't provide this derivative
%                function, please input "Null".
%   niter :      the maximum number of iterations
%   par:         a cell with two elements, the first denotes whether we need 
%                all the estimates in whole process of just the final estimate 
%                we obtain, "1" means need, "0" means not; the second
%                denotes how many iteration you want AMP to run, if you
%                input a positive integer t, then AMP runs t times
%                iterations for you, if you input the string 'Auto', then
%                AMP will try to runs 100 times iterations and stop when
%                the ratio of l2 norm of (x_(t+1)-x_(t)) and l2 norm of
%                x_(t) less then 0.01
% Related functions: Eta_der_Estimate_CAMP

% check: @A,@Eta,@Etader, sigma is the estimated std instead of variance.
% sigma_w seems useless in this function

warning off;
t1     = tic;
n = length(y);
lengthN = A'*zeros(n,1);
N = length(lengthN);

par1=logical(par{1});
par2=par{2};
if ischar(par2)
    niter=100;
else
    niter=par2;
end

empiricaliterwatch_sigma=zeros(niter+1,1);
xall=zeros(N,niter+1);
mx=zeros(N,1);
mz=y-A*(mx./colnormA);
% iteration_time_total=0;

t2     = tic;
for iter=1:niter
    temp_z=A'*(mz)./colnormA+mx;
    sigma_hat= 1/sqrt(log(2))*median(abs(temp_z));
    mx=Eta(temp_z,sigma_hat);
    
    if strcmpi(Etader,'Null')
        mz=y-A*(mx./colnormA)+mz*Eta_der_Estimate_CAMP(temp_z,sigma_hat,Eta )*N/n;
    else
        [etaderR,etaderI]=Etader(temp_z,sigma_hat);
        mz=y-A*(mx./colnormA)+mz*(sum(etaderR)+sum(etaderI))/(2*n);
    end
    
    empiricaliterwatch_sigma(iter)=sigma_hat;
    xall(:,iter+1)=mx./colnormA;
    tn=mx./colnormA;
    Axb=y-A*tn;
    
    obj=sum(conj(Axb).*Axb)/2;
    out1.obj = obj;
    Obj(iter)= obj;
    time2        = toc(t2);
    %      error = abs( empiricaliterwatch_sigma(iter)- ...
    %         1/sqrt(log(2))*median(abs(A'*(mz)./colnormA+mx))); 
    if niter==100 && abs( empiricaliterwatch_sigma(iter)- ...
            1/sqrt(log(2))*median(abs(A'*(mz)./colnormA+mx)) )<0.0001
        break;
    end
    
end
time        = toc(t1);

out1.time    = time;
out1.time2   = time2;
empiricaliterwatch_sigma(iter+1)=1/sqrt(log(2))*median(abs(A'*(mz)./colnormA+mx));

if niter==100
    empiricaliterwatch_sigma = empiricaliterwatch_sigma(1:(iter+1));
    xall=xall(:,1:(iter+1));
end

if iter==100
    %  fprintf('Iteration reaches the maximum (100) times,\nthe algorithm does not converge within 100 iterations.\n');
end

% if par1
%     xhat=xall;
% else
%     xhat=xall(:,end);
% end

out1.sol  = xall(:,end);
out1.iter = iter;
out1.obj  = Obj(iter);
end