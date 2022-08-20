clc; clear; close all;
path(path,'/data/home/u19121573/afternorm/demonHW_xujing');
n           = 2048;
m           = 839; 
s           = 20;  

start = tic;
fprintf(' Start to generate the compressed sensing data...\n'); 
 
%%  Huawei A*(The correlation is 0.0345)  
A   = zeros(m,n);
for m1 = 1:839
    for n1 = 1:832
        A(m1,n1) = exp(1i*420*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n1-1)/839);
    end
    for n2 = 833:1664
        A(m1,n2) = exp(1i*419*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n2-833)/839);
    end
    for n3 = 1665:2048
        A(m1,n3) = exp(1i*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n3-1665)/839);
    end
end
%% Huawei A**(839x5952，The correlation is 0.0345 )    
% A   = zeros(m,n);
% for m1 = 1:839
%     for n1 = 1:837
%         A(m1,n1) = exp(1i*420*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n1-1)/839);
%     end
%     for n2 = 838:1674
%          A(m1,n2) = exp(1i*419*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n2-838)/839);
%     end
%     for n3 = 1675:2511
%         A(m1,n3) = exp(1i*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n3-1675)/839);
%     end
%     for n4 = 2512:3348
%         A(m1,n4) = exp(1i*838*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n4-2512)/839);
%     end
%     for n5 = 3349:4185
%         A(m1,n5) = exp(1i*15*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n5-3349)/839);
%     end
%     for n6 = 4186:5022
%          A(m1,n6) = exp(1i*824*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n6-4186)/839);
%     end
%     for n7 = 5023:5859
%         A(m1,n7) = exp(1i*427*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n7-5023)/839);
%     end
%     for n8 = 5860:5952
%         A(m1,n8) = exp(1i*412*pi*m1*(m1-1)/839)*exp(1i*2*pi*(m1-1)*(n8-5860)/839);
%     end
% end

%% 
A = A./sqrt(839);
data.A  = A;       % 无normalization步骤

% n00 = randperm(64,s);  % 产生固定激活用户
load('n00', 'n00');
Reaj = zeros(1,64);    % 初始化均为非激活状态
Reaj(n00) = 1;         % 其中1表示用户激活，0表示用户未激活

data.At = data.A';                                
clear Ax Ay m1 n1 m2 n2
fprintf(' Data generation used %2.4f seconds.\n\n',toc(start)); 

pars.draw = 1;
pars.eta  = 1;

% 重复仿真计算虚警概率概率
Iter  = 200 ;
OBJ   = zeros(Iter,1);
iter  = zeros(Iter,1);
Time  = zeros(Iter,1);
Err   = zeros(Iter,1);

b  = zeros(m,1) ;

fprintf(' \n 组别   方差    CPU时间     相对误差     目标值    迭代次数    虚警概率    漏检概率  \n');
fprintf('\n ----------------------------------------------------------------------------------\n');
Var =[0.005,0.0075,0.02,0.06];
% Var =[0.000,0.001,0.005,0.010,0.015,0.02,0.025,0.030,0.040,0.045,0.050,0.060,0.070,0.080,0.090,0.100];
% Var =[0.000,0.001,0.005,0.010,0.015,0.02,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,,0.090,0.100];
% Var =[0.01,0.05,0.10,...
%     0.15, 0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,...
%     0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,...
%     1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.25,...
%     3.50,3.75,4.00,4.25,4.50,4.75,5.00];
Threshold(1)  = 0.21;
Threshold(2)  = 0.22;
Threshold(3)  = 0.41;
Threshold(4)  = 0.72;
xbeta  = sqrt(0.5*1);
for v=1:length(Var)
    var = Var(v);
    RA_C = zeros(1,64);
    for i =1: Iter
        x0     = zeros(32,64);
        for n0 = n00
            m0 = randperm(32,1);
            x0(m0,n0) = xbeta*randn+1i*xbeta*randn;
        end
        x0     = x0(:);
     %   data.x_opt = x0;   
        noise_sigma = sqrt(0.5*var);
        noise   = noise_sigma*randn(m,1)+ 1i*noise_sigma*randn(m,1);
        b  = A*x0 + noise; 
        data.b  = b;               %  修改部分
        out1   = BNHTP(n,s,data,pars);
       % x00  = out1.sol;  
        x00   = out1.sol.*(abs(out1.sol)>Threshold(v));   % 阈值处理
        Rec   = reshape(x00,32,64);       % 恢复值
        Recj  = (sign(sum(abs(Rec)))) ;   % 恢复值的激活用户集 
        Ra_c  = Reaj-Recj  ;              % 如果元素是+1就是漏检，如果是-1就是虚警     
        RA_C  = RA_C+Ra_c; 
        OBJ(i)  = out1.obj;
        iter(i) = out1.iter;
        Time(i) = out1.time;
        Err(i)  = norm(x0-x00)/norm(x0);
    end
                 
rA_C  = RA_C/Iter ;                      
Fal   = abs( sum(rA_C(rA_C<0))/(64-s)) ; 
Mis   = sum(rA_C(rA_C>0))/(s) ;          

fprintf('\n   %d',  v);
fprintf('   %.3f % ', Var(v));
fprintf('   %.3f ',  sum(Time)/Iter);
fprintf('     %5.2e',...
    sum(Err)/Iter);
fprintf('    %5.2e', sum(OBJ)/Iter);
fprintf('    %.2f',  sum(iter)/Iter);
fprintf('      %.4f % ',  Fal);
fprintf('    %.4f % \n', Mis);
end
fprintf('\n ----------------------------------------------------------------------------------\n');

TIME = toc(start);

