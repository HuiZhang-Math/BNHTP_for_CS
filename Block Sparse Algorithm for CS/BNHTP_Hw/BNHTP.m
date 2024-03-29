function Out = BNHTP(n, s, data, pars)

%     data    : A structure (data.A, data.At, data.b), where       (required)
%               data.A \in C^{m-x-n}   
%               data.At is the conjugate transpose of data.A 
%               data.b \in C^{m-x-1} 
%     pars:     Parameters are all OPTIONAL(required)
%               pars.x0      --  Starting point of x,   pars.x0=zeros(n,1) (default)
%               pars.eta     --  A positive parameter,  a default one is given related to inputs  
%               pars.display --  =1. Display results for each iteration.(default)
%                                =0. Don't display results for each iteration.
%               pars.draw    --  A  graph will be drawn if pars.draw=1 (default) 
%                                No graph will be drawn if pars.draw=0
%               pars.maxit   --  Maximum number of iterations, (default,2000) 
%               pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%
% Outputs:
%     Out.sol:           The sparse solution x
%     Out.sparsity:      Sparsity level of Out.sol
%     Out.normgrad:      L2 norm of the gradient at Out.sol  
%     Out.error:         Error used to terminate this solver 
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%     Out.grad:          Gradient at Out.sol
%     Out.obj:           Objective function value at Out.sol 
%
% This code is programmed based on the algorithm proposed in 
% "S. Zhou, N. Xiu and H. Qi, Global and Quadratic Convergence of Newton 
% Hard-Thresholding Pursuit, Journal of Machine Learning Research, 2019."


warning off;
t0     = tic;

if     nargin<2
       disp(' No enough inputs. No problems will be solverd!'); return;
end

 func      = @(x,fgh,T1,T2)csfun(x,fgh,T1,T2,data); 
if nargin>=2
    if nargin<3; pars=[]; end
    if isfield(pars,'display');display = pars.display;else; display = 1;    end
    if isfield(pars,'draw');   draw    = pars.draw;   else; draw    = 1;    end
    if isfield(pars,'maxit');  itmax   = pars.maxit;  else; itmax   = 2000; end
    if isfield(pars,'tol');    tol     = pars.tol;    else; tol     = 1e-6; end  
    if isfield(pars,'x0');     x0      = pars.x0;     else; x0 = zeros(n,1);end 
    
   [x0,obj,g,eta] = getparameters(n,s,x0,func,pars);     
end
   
x       = x0;
beta    = 0.5;
sigma   = 5e-5;
I       = 1:n;
delta   = 1e-10;
pcgtol  = 0.01*tol*s;
T0      = zeros(s,1);
Error   = zeros(1,itmax);
Obj     = zeros(1,itmax);
dg1     = zeros(1,itmax);
FNorm   = @(x)norm(x)^2;

% if display 
% fprintf(' Start to run the solver...\n'); 
% fprintf('\n Iter             Error                Ojective \n'); 
% fprintf('------------------------------------------------\n');
% end

% Initial check for the starting point
if FNorm(g)==0 && nnz(x)<=s
fprintf('Starting point is a good stationary point, try another one...\n'); 
Out.sol = x;
Out.obj = obj;
return;
end

if max(isnan(g))
x0      = zeros(n,1);
rind    = randi(n);
x0(rind)= rand;
[obj,g] = func(x0,'ObjGrad',[],[]);
end
 
d0      = -g;       
% The main body  
for iter = 1:itmax
    
%% �ֿ���֧�ż�
    TI    = [];
    xtg   = x0-eta*g;
    for  i = 1:64
        [~,Ti] = maxk(abs(xtg((32*(i-1)+1):(32*i))),1);
        TI     = union(TI, Ti + 32 * (i-1));
    end
    TIc    = setdiff(I,TI);
    Xtg    = abs(xtg);
    Xtg(TIc) = 0;
    [~,T] = maxk(abs(Xtg),s); 
    
    flag  = isempty(setdiff(T,T0));          
    Tc    = setdiff(I,T);
    
    % Calculate the error for stopping criteria   
    xtaus       = max(0,max(abs(g(Tc)))-min(abs(x(T)))/eta);
    if flag
    FxT         = sqrt(FNorm(g(T)));
    Error(iter) = xtaus + FxT;
    else
    FxT         = sqrt(FNorm(g(T))+ FNorm(x(Tc)));
    Error(iter) = xtaus + FxT;    
    end
     
%     if display 
%     fprintf('%4d             %5.2e              %5.2e\n',iter,Error(iter),obj); 
%     end
             
    % Stopping criteria
    if Error(iter)<tol; break;  end
   
    % update next iterate
    if  iter   == 1 || flag           % update next iterate if T==supp(x^k)
        H       =  func(x0,'Hess',T,[]);
        if ~isa(H,'function_handle')
            d   = -H\g(T);
        else
            [d,~]= pcg(H,-g(T),pcgtol,20);
        end
        dg      = sum(d.*g(T));
        ngT     = FNorm(g(T));
        if dg   > max(-delta*FNorm(d), -ngT) || isnan(dg)
            d = -g(T);
            d = d + 1/32 * d0(T);
            d0(T) = d;
            dg  = ngT;
        end
    else                              % update next iterate if T~=supp(x^k)
        TTc     = intersect(T0,Tc); 
        [H,D]   = func(x0,'Hess',T,TTc);
        
        if isa(D,'function_handle')
           Dx   = D(x0(TTc));
        else
           Dx   = D*x0(TTc);
        end
        
        if ~isa(H,'function_handle')
            d   = H\( Dx- g(T));
        else
           [d,~]= pcg(H,Dx- g(T), pcgtol,20);           
        end
        
        dgT     = sum(d.*g(T));
        dg      = dgT-sum(x0(TTc).*g(TTc));
%         Fnz     = FNorm(x(TTc))/4/eta;
%         delta0  = delta;
%         if Fnz  > 1e-4; delta0 = 1e-4; end
%         ngT     = FNorm(g(T));
%         if dgT  > max(-delta0*FNorm(d)+Fnz, -ngT) || isnan(dg) 
%         d       = -g(T); 
%         dg      = ngT; 
%         end            
    end
    
    alpha    = 1; 
    x(Tc)    = 0;    
    obj0     = obj;        
    Obj(iter)= obj;
    dg1(iter)= dg;
     
    % Amijio line search
    for i      = 1:6
        x(T)   = x0(T) + alpha*d;
        obj    = func(x,'ObjGrad',[],[]);
        if obj < obj0  + alpha*sigma*dg; break; end        
        alpha  = beta*alpha;
    end
    
    % Hard Thresholding Pursuit if the obj increases
    fhtp    = 0;
    if obj  > obj0 
       x(T) = xtg(T); 
       obj  = func(x,'ObjGrad',[],[]); 
       fhtp = 1;
    end
    
    % Stopping criteria
    flag1   = (abs(obj-obj0)<1e-3*(1+abs(obj)) && fhtp);               % pre 1e-6
    flag2   = (abs(obj-obj0)<1e-3*(1+abs(obj))&& Error(iter)<1e-2);    % pre 1e-10
    if  iter>10 &&  (flag1 || flag2)      
        if obj > obj0
           iter    = iter-1; 
           x       = x0; 
           T       = T0; 
        end   
        break;
     end 
 
    T0      = T; 
    x0      = x; 
    [obj,g] = func(x,'ObjGrad',[],[]);
    
    
    % Update eta
    if  mod(iter,20)==0  
        if Error(iter)>1/iter^2  
        if iter<1500; eta = eta/1.25; 
        else;         eta = eta/1.5; 
        end     
        else;         eta = eta*1.25;   
        end
    end     
end

[obj ,g]    = func(x,'ObjGrad',[],[]);
% results output
time        = toc(t0);
Out.normgrad= sqrt(FNorm(g)); 
Out.error   = sqrt(FNorm(g(T))+ FNorm(x(setdiff(I,T))));
Out.time    = time;
Out.iter    = iter;
Out.grad    = g;
Out.sol     = x;
Out.obj     = obj; 
Out.Error   = Error;
Obj(iter)   = obj; 
Out.Obj     = Obj;
Out.dg1     = dg1;

% if draw
%     figure
%     subplot(121) 
%     Obj(iter)= obj;  
%     PlotFun(Obj,iter,'r.-','f(x^k)'); 
%     subplot(122) 
%     PlotFun(Error,iter,'r.-','error') 
% end


if display 
  % fprintf('------------------------------------------------\n');
   if Out.normgrad<1e-5
     % fprintf(' A global optimal solution might be found\n');
    %  fprintf(' because of ||gradient||=%5.2e!\n',Out.normgrad); 
      if Out.iter>1500
      fprintf('\n Since the number of iterations reaches to %d\n',Out.iter);
      fprintf(' Try to rerun the solver with a smaller pars.eta=%5.2e\n',Out.eta); 
      end
   %   fprintf('------------------------------------------------\n');
   end
end


end

%--------------------------------------------------------------------------
function [x0,obj,g,eta]=getparameters(n,s,x0,func,pars)

    if isfield(pars,'x0') && norm(x0)>0
       [obj0,g0] = func(zeros(n,1),'ObjGrad',[],[]);  
       [obj,g]   = func(pars.x0,'ObjGrad',[],[]); 
       if obj0   < obj/10
          x0     =  zeros(n,1); 
          obj    = obj0;  
          g      = g0; 
       else  
          ns0    = nnz(pars.x0);
          if ns0==s
          [~,T]    = maxk(pars.x0,s,'ComparisonMethod','abs'); 
          x0       = pars.x0;  
          pars.eta = min(abs(x0(T)))/(1+max(abs(g(setdiff(1:n, T)))));   
          elseif ns0<s
          x0        = pars.x0;  
          pars.eta  = max(x0(x0>0.1))/(1+max(abs(g)));   
          else 
          [~,T]     = maxk(pars.x0,s,'ComparisonMethod','abs'); 
          x0        = zeros(n,1);
          x0(T)     = pars.x0(T);  
          pars.eta  = max(x0(x0>0.1))/(1+max(abs(g)));  
          end
          
          if isempty(pars.eta) 
          pars.eta  = max(abs(x0))/(1+max(abs(g))); 
          end
          
       end
    else
        [obj,g]  = func(x0,'ObjGrad',[],[]); 
    end
 
    
    if isfield(pars,'eta')      
        eta  = pars.eta;       
    else % set a proper parameter eta
        [~,g1] = func(ones(n,1),'ObjGrad',[],[]) ;
        abg1   = abs(g1);
        T      = find(abg1>1e-8);
        maxe   = sum(1./(abg1(T)+eps))/nnz(T);
        if  isempty(T) 
            eta    = 10*(1+s/n)/min(10, log(n));
        else
            if maxe>2
            eta  = (log2(1+maxe)/log2(maxe))*exp((s/n)^(1/3));
            elseif maxe<1
            eta  = (log2(1+ maxe))*(n/s)^(1/2);    
            else
            eta  = (log2(1+ maxe))*exp((s/n)^(1/3));
            end     
        end
    end 
    
end
%--------------------------------------------------------------------------
function [out1,out2] = csfun(x,fgh,T1,T2,data)
    Tx = find(x); 
    Ax = data.A(:,Tx)*x(Tx);
    switch fgh        
    case 'ObjGrad'
        Axb  = Ax-data.b;
        out1 = sum(conj(Axb).*Axb)/2;         %objective function value of f
        if  nargout>1 
        out2 = data.At*Axb;                    %gradien of f
        end
    case 'Hess'
        
        if  length(T1)<2000
            out1 = data.At(T1,:)*data.A(:,T1);       %submatrix containing T1 rows and T1 columns of Hessian
        else
            AT1  = data.A(:,T1);
            out1 = @(v)((AT1*v)'*AT1)';     
        end
        
        if nargout>1
        out2 = @(v)(data.At(T1,:)*(data.A(:,T2)*v));  %submatrix containing T1 rows and T2 columns of Hessian
        end  
        
    end
end

%--------------------------------------------------------------------------
% function  PlotFun(input,iter,c, key) 
%     if  input(iter)>1e-40 && input(iter)<1e-5
%         semilogy(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
%     else
%         plot(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
%     end
%     xlabel('Iter'); ylabel(key); grid on    
%     
% end

