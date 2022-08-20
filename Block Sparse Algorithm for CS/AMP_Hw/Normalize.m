function [A,colnormA] = Normalize(y,A)

n = length(y);
lengthN = A'*zeros(n,1);
N = length(lengthN);
%%%%% Normalize A
disp('Normalizing A matrix, please wait...');

pick=randperm(N);
DtmIndex=pick(1:5);
DtmNorms=zeros(5,1);
I=eye(N);

for i=1:5
    % DtmNorms(i)=sum(abs( A(I(:,DtmIndex(i)),1) ).^2);
    DtmNorms(i)=sum(abs( A*(I(:,DtmIndex(i)))).^2);
end

if (sum(DtmNorms>1.1)+sum(DtmNorms<0.9))>0 % We need to normalize A matrix
   % disp('It is necessary to normalize the A matrix, please wait...');
    tempA=zeros(n,N);
    tempA_ra=zeros(n,N);
    normalize_time_total=0;
    for j=1:N
       t0=cputime;
       %  tempA(:,j)=A(I(:,j),1);  
        tempA(:,j)=A*(I(:,j));  
        normalize_time=(cputime-t0)/60;
        normalize_time_total=normalize_time_total+normalize_time;
        
        if j/(N/100)==fix(j/(N/100))
            normalize_time_remain=normalize_time_total*(N-j)/j;
            percent=j/(N/100);
     %       disp(['Normalizing has been through ' num2str(percent) '%.' 10 'The estimated remaining time for Normalizing is ' num2str(normalize_time_remain) ' minutes.']);
        end
        
    end

    for j=1:N                              %remove average
        tempA_ra(:,j)=tempA(:,j)-mean(tempA(:,j));
    end
     colnormA=(sqrt(sum(abs(tempA_ra).^2,1)))';
    ind=find(colnormA==0);
    colnormA(ind)=(sqrt(sum(abs(tempA(:,ind)).^2,1)))';
    
  %  disp('Normalizing ends, Iteration starting...');
else
    disp('It is not necessary to normalize the A matrix, Iteration starting...');
    colnormA=ones(N,1);
end
end


