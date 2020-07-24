function [P,E,B]=LAM_learning(F0,W0,lambda,alpha)
% This function inputs the initial PCM (F0) and an affinity matrx (W0),
% and outputs an estimated LAM (P).


tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;

n=length(F0);


B=zeros(n);

E=B;
Lambda1=B;
Lambda2=B;
Lambda3=B;
C=B;

L=diag(sum(W0))-W0;

L=sparse(L);
I=eye(n);



for iter=1:max_iter
%     iter
    % for the B subproblem (actually the )
    Temp1=W0-E+Lambda1/mu;
    Temp2=B-Lambda2/mu;
    Temp3=C-Lambda3/mu;
    Temp=(Temp1+Temp2+Temp3)/3;
    Temp=(Temp+Temp')/2;
    P=prox_nuclear(Temp,1/(mu*3));
    
    
    % for the E subproblem
%     Temp=W0-W+Lambda1/mu;
%     E=prox_l1(Temp,(1*lambda)/(mu));
    Temp=W0-P+Lambda1/mu;
    E=prox_l1(Temp,(1*lambda)/(mu));
    
   
    
    % The B subproblem
    Temp=P+Lambda2./mu;
   
    B=Temp;
%     D=F-Lambda3./mu;
    B(F0==1)=1;
    B(F0==-1)=0;
    B(B<0)=0;
    B(B>1)=1;
    
    % The C subproblem
    C=(mu*P+Lambda3)*inv(2*alpha*L+mu*I);    
   
    
    d1=W0-P-E;
    d2=(P-B);
    d3=P-C;
    
    chg = max([ max(abs(d1(:))),max(abs(d2(:))),max(abs(d3(:))) ]);
    if chg < tol
        break;
    end 
    
    Lambda1=Lambda1+mu*d1;
    Lambda2=Lambda2+mu*d2;
    Lambda3=Lambda3+mu*d3;
    mu = min(rho*mu,max_mu); 
end