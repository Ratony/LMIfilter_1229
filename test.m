%function [Af,Bf,Cf,Lam,Gamma,W,flag]=LMI_FiniteTime_Filter(A,B,C,D,L,pro,pro2,taoM,dM,kesai,Ptrans)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper

clear all
clc
% Input:

% Output:
%known parameters
k1=0.1;
k2=0.1;
delta=0.2;
A=cell(1,2);
B=cell(1,2);
C=cell(1,2);
L=cell(1,2);

A{1}=[-1 0 0;0 -2 1;0 0 -1];
B{1}=[1 0 0.3]';
C{1}=[0 0.1 0.3];
L{1}=[0.5 -0.1 0.1];

A{2}=[-2 0.3 0;0 -2 0;0 0 -1.8];
B{2}=[-0.6 0.5 0]';
C{2}=[0.2 0.2 0];
L{2}=[0 -0.3 0.2];

% if 1
%     A{1}=[-9.1 50;-1 -10];
%     B{1}=[0 1]';
%     C{1}=[1 0];
%     L{1}=[1 0];
%     
%     A{2}=[-0.1 50;-1 -16];
%     B{2}=[0 1]';
%     C{2}=[1 0];
%     L{2}=[1 0];
% else
%     
%     A{1}=[-1 2;-3 0.1];
%     B{1}=[-0.2 1]';
%     C{1}=[1 0.5];
%     L{1}=[0.6 1];
%     
%     A{2}=[-1 2;-3 0.1];
%     B{2}=[-0.2 1]';
%     C{2}=[1 0.5];
%     L{2}=[0.6 1];
% end
pro1=0.2;
pro2=0.8;
beta=0.8;
%% Variables assignment
[m,n]=size(B{1})
r=n+1;
m1=m*2;
R=1*eye(m1);

%% constant value
T1=0.5;
I1=0.3;
T2=1/T1;
I2=1/I1;
%% SDP variables

s1=sdpvar(m1);
s2=sdpvar(m1);

w1=sdpvar(m1);
w2=sdpvar(m1);

P1=sdpvar(m);
P2=sdpvar(m);
% P1=cell(1,2);
% for i=1:2
%     P1{i}=sdpvar(m);
% end

% for i=1:4
%     Q{i}=sdpvar(m);
% end


%% performance index
Gamma=sdpvar(n);

%% event-triggered scheme parameters
Lam1=sdpvar(1);
%Lam1=sdpvar(1);


%% filter variables and P
Af=cell(2,2);
Bf=cell(2,2);
Cf=cell(2,2);
for i=1:2
    for j=1:2
            Af{i}{j}=sdpvar(m);
            Bf{i}{j}=sdpvar(m,n);
            Cf{i}{j}=sdpvar(n,m);
    end
end

% s1=R^(1/2)*S1*R^(1/2);
% s2=R^(1/2)*S2*R^(1/2);
% w1=R^(1/2)*W1*R^(1/2);
% w2=R^(1/2)*W2*R^(1/2);
Q=[P1 -P2;-P2 P2];
%[-pro1* e1;pro1*Bf{1}*C{1} e1];
e1=zeros(m);
%% Theorem 2 LMIs
flag=0;
Phi1=cell(2,2);
for i=1:2
    for j=1:2
            %%%%%i表示对象的下标取值，j表示滤波器下标取值；
            %定义
            Phi1{i}{j}=blkvar;
            
            PA=[P1*A{i} -Af{i}{j};-P2*A{i} Af{i}{j}];
            PB1=[-pro1*Bf{i}{j}*C{i} e1;pro1*Bf{i}{j}*C{i} e1];   
            PB2=[-pro2*Bf{i}{j}*C{i} e1;pro2*Bf{i}{j}*C{i} e1];
            PB3=[-pro2*Bf{i}{j};pro2*Bf{i}{j}];
            PB4=[P1*B{i};-P2*B{i}];
            Eij=[L{i}';-Cf{i}{j}'];
            Cij=[delta*C{i}'*Lam1*C{i} e1;e1 e1];
            
            
            %列举
            Phi1{i}{j}(1,1)=PA+PA'+s1+s2-beta*Q-T2*w1-I2*w2;
        
            Phi1{i}{j}(1,2)=PB1+T2*w1; 
            Phi1{i}{j}(1,4)=I2*w2+PB2;
            Phi1{i}{j}(1,6)=-PB3;
            Phi1{i}{j}(1,7)=PB4;
            
            Phi1{i}{j}(1,8)=T1*PA';
            Phi1{i}{j}(1,9)=I1*PA';
            Phi1{i}{j}(1,10)=Eij;
      
            Phi1{i}{j}(2,2)=-exp(beta*T1)*(1-k1)*s1-2*T2*w1;
            Phi1{i}{j}(2,3)=T2*w1;
            Phi1{i}{j}(2,8)=T1*PB1';
            Phi1{i}{j}(2,9)=I1*PB1';
         
            
            Phi1{i}{j}(3,3)=-T2*w1;
                       
            
            Phi1{i}{j}(4,4)=-exp(beta*I1)*(1-k2)*s2-2*T2*w2+Cij;
            Phi1{i}{j}(4,8)=T1*PB2';
            Phi1{i}{j}(4,9)=I1*PB2';
            
            Phi1{i}{j}(4,5)=I2*w2;
            
            Phi1{i}{j}(5,5)=-I2*w2;
       
            Phi1{i}{j}(6,6)=-Lam1;
            Phi1{i}{j}(6,8)=-PB3';
            Phi1{i}{j}(6,9)=-PB3';
            
            Phi1{i}{j}(7,7)=-Gamma;
            
            Phi1{i}{j}(7,8)=T1*PB4';
            Phi1{i}{j}(7,9)=I1*PB4';
         
            Phi1{i}{j}(8,8)=T1*(w1-Q-Q');
        
            Phi1{i}{j}(9,9)=I1*(w2-Q-Q');  
            Phi1{i}{j}(10,10)=-eye(1);     
            
            Phi1{i}{j}=sdpvar(Phi1{i}{j});
    end
end



%constrain=max(eig(Q1))
%% solution of LMIs

LMIs=[Q>0,s1>0,s2>0,w1>0,w2>0,1>Gamma>0,Lam1>0,Phi1{1}{1}<0,...,
    Phi1{1}{2}<0,Phi1{2}{1}<0,Phi1{2}{2}<0];
%sol=solvesdp(LMIs,Gamma)
if 1
options=sdpsettings('solver','sdpt3','verbose',0);
%sol=solvesdp(LMIs)
sol=optimize(LMIs)
sol.solvertime
sol.yalmiptime
end
if 0
options=sdpsettings('solver','sdpt3','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0
if sol.problem == 0 
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end
end
%% optimize
if 0
    options=sdpsettings('solver','sedumi','verbose',0);
    sol=optimize(LMIs,[],options); 

    if sol.problem == 0
            valobj=value(Phi1);
            if 0
                Gamma=check(LMIs);
                if min(Gamma)<1
                    Gamma=Gamma;
                end
            end
    else
        disp('error');
    end

end

%% solve unknown decesion variables
if 1
    
for i=1:2
    for j=1:2
        Af{i}{j}=value(Af{i}{j});
        Bf{i}{j}=value(Bf{i}{j});
        Cf{i}{j}=value(Cf{i}{j});
    end
end

P1=value(P1)
P2=value(P2)
Lam1=value(Lam1);
%Lam1=value(Lam1);
Gamma=value(Gamma);
%xlswrite('filter_parameter.xls',Af);
end
%end