%function [ET_TIME,ET_INTER]=Pro_polot_event(A,B,C,D,L,AF,BF,CF,Lam,KESAI1,X0,XF0,Exp_Dimensions)
%h1=0.002;
%h2=0.35;

%KESAI1=0.2;

%%3
%X0=[1 1 -2]';
%XF0=[0 0 0]';

%XF0=[0.1 0.05 0.1]';
clear all
clc
[m,AF,BF,CF,delta,A,B,C,L,Lam1,R]=test()
Exp_Dimensions=m;
W0=0;
    belta1l=0.5;
    belta2l=0.5;
    belta1u=0.5;
    belta2u=0.5;
D=-1;   
X0=[-0 0 0]';
XF0=[0 0 0]';


    w1g(1:1,1)=(1/(exp(4)));
    w2g(1:1,1)=1-w1g(1:1,1);
    
    m1g(1:1,1)=0.2;
    m2g(1:1,1)=0.8;

  
     
h=0.01;
h=0.01;

finaltime=2;


W0=0;
DX0=w1g(1:1,1)*A{1}*X0+w1g(1:1,1)*B{1}*W0+w2g(1:1,1)*A{2}*X0+w2g(1:1,1)*B{2}*W0;
X(:,1)=DX0*h+X0;
Z0=w1g(1:1,1)*L{1}*X0+w2g(1:1,1)*L{2}*X0;
Y0=w1g(1:1,1)*C{1}*X(:,1)+w2g(1:1,1)*C{2}*X(:,1)+D*W0


DXF0=m1g(1:1,1)*AF{1}*XF0+m1g(1:1,1)*BF{1}*Y0+m2g(1:1,1)*AF{2}*XF0+m2g(1:1,1)*BF{2}*Y0;
XF(:,1)=DXF0*h+XF0;
ZF0=m1g(1:1,1)*CF{1}*XF0+m2g(1:1,1)*CF{2}*XF0;


%��һ��
ET_TIME(1,1)=0;
%��������
TIMES=1;
%����ֵ
ET_Value(:,TIMES)=Y0; 
%��������
ET_INTER(1,1)=0;

Time_tri=0;
Loss=0;
Bernoulli_Event=binornd(1,0.2,1,finaltime/h);

for COUNT=1:finaltime/h
    
   
    if COUNT==1
        Bernu(:,COUNT)=0;
    else
        Bernu(:,COUNT)=Bernoulli_Event(COUNT);
    end

    t=COUNT*h;
%% �˲�����ʱ�������Ⱥ���
     m1g(:,COUNT)=exp(-X(1:1,COUNT)^2/8);
     m2g(:,COUNT)=1-m1g(:,COUNT);
   
%%  �����ʱ�������Ⱥ���   
    w1g(:,COUNT)=sin(X(1:1,COUNT))^2;
    w2g(:,COUNT)=1-w1g(1:1,COUNT);
    m1g(:,COUNT)=cos(X(1:1,COUNT))^2;
    m2g(:,COUNT)=1-m1g(:,COUNT);
    
    W(:,COUNT)=sin(4*pi*t)*exp(-2*t);
    
  
    DX(:,COUNT)=w1g(:,COUNT)*A{1}*X(:,COUNT)+w1g(:,COUNT)*B{1}*W(:,COUNT)+w2g(:,COUNT)*A{2}*X(:,COUNT)+w2g(:,COUNT)*B{2}*W(:,COUNT);
    X(:,COUNT+1)=DX(:,COUNT)*h+X(:,COUNT);
    Z(:,COUNT)=w1g(:,COUNT)*L{1}*X(:,COUNT)+w2g(:,COUNT)*L{2}*X(:,COUNT);
    Y(:,COUNT)=w1g(:,COUNT)*C{1}*X(:,COUNT+1)+w2g(:,COUNT)*C{2}*X(:,COUNT+1)+D*W(:,COUNT);
    
    if Bernu(:,COUNT)==0
         if((Y(:,COUNT)-ET_Value(:,TIMES))'*Lam1*(Y(:,COUNT)-ET_Value(:,TIMES)))>(delta*ET_Value(:,TIMES)'*Lam1*ET_Value(:,TIMES))
            TIMES=TIMES+1;
            ET_Value(:,TIMES)=Y(:,COUNT);
            ET_TIME(TIMES,1)=t;
            ET_INTER(TIMES,1)=ET_TIME(TIMES,1)-ET_TIME(TIMES-1,1);

            DXF(:,COUNT)=m1g(:,COUNT)*AF{1}*XF(:,COUNT)+m1g(:,COUNT)*BF{1}*ET_Value(:,TIMES)+m2g(:,COUNT)*AF{2}*XF(:,COUNT)+m2g(:,COUNT)*BF{2}*ET_Value(:,TIMES);
            XF(:,COUNT+1)=DXF(:,COUNT)*h+XF(:,COUNT);
            ZF(:,COUNT)=m1g(:,COUNT)*CF{1}*XF(:,COUNT)+m2g(:,COUNT)*CF{2}*XF(:,COUNT);
            disp('test')
        else
            DXF(:,COUNT)=m1g(:,COUNT)*AF{1}*XF(:,COUNT)+m1g(:,COUNT)*BF{1}*ET_Value(:,TIMES)+m2g(:,COUNT)*AF{2}*XF(:,COUNT)+m2g(:,COUNT)*BF{2}*ET_Value(:,TIMES);
            XF(:,COUNT+1)=DXF(:,COUNT)*h+XF(:,COUNT);
            ZF(:,COUNT)=m1g(:,COUNT)*CF{1}*XF(:,COUNT)+m2g(:,COUNT)*CF{2}*XF(:,COUNT);
        end
    else
        Time_tri= Time_tri+1    ;     
        DXF(:,COUNT)=m1g(:,COUNT)*AF{1}*XF(:,COUNT)+m1g(:,COUNT)*BF{1}*ET_Value(:,TIMES)+m2g(:,COUNT)*AF{2}*XF(:,COUNT)+m2g(:,COUNT)*BF{2}*ET_Value(:,TIMES);
        XF(:,COUNT+1)=DXF(:,COUNT)*h+XF(:,COUNT);
        ZF(:,COUNT)=m1g(:,COUNT)*CF{1}*XF(:,COUNT)+m2g(:,COUNT)*CF{2}*XF(:,COUNT);
    end
   E(:,COUNT)=Z(:,COUNT)-ZF(:,COUNT);
   EXT(:,COUNT)=[X(:,COUNT);XF(:,COUNT)];
   FINITE(:,COUNT)=EXT(:,COUNT)'*R*EXT(:,COUNT);
   c(:,COUNT)=0.5;
end

A=TIMES
Time_tri
figure(1);
plot([0:h:finaltime-h],E(:),'g-',[0:h:finaltime-h],Z(:),'r-',[0:h:finaltime-h],ZF(:),'k-')
legend('e(t)','z(t)','zf(t)');
grid on;
figure(2);
stem(ET_TIME(:,1),ET_INTER(:,1));
 
if Exp_Dimensions==3 
figure(3);
plot([0:h:finaltime],X(1,:),'k',[0:h:finaltime],XF(1,:),'k:',[0:h:finaltime],X(2,:),'r',[0:h:finaltime],XF(2,:),'r:',[0:h:finaltime],X(3,:),'g',[0:h:finaltime],XF(3,:),'g')
legend('x1(t)','xf1(t)','x2(t)','xf2(t)','x3(t)','xf3(t)');
else
figure(3);
plot([0:h:finaltime],X(1,:),'k',[0:h:finaltime],XF(1,:),'k:',[0:h:finaltime],X(2,:),'r',[0:h:finaltime],XF(2,:),'r:')
legend('x1(t)','xf1(t)','x2(t)','xf2(t)');
end

figure(4);
plot([0:h:finaltime-h],FINITE(1,:),'k',[0:h:finaltime-h],c(:,COUNT),'m')
legend('bound index')

%xlable('x(t))');
if 0
    figure(4);
    plot([0:h:finaltime],X(2,:),'k-',[0:h:finaltime],XF(2,:),'k.')
    legend('x(2,:)','xf(2,:)')
    grid on; 
else
    figure(5);
    plot([0:h:finaltime],X(1,:),'k',[0:h:finaltime],X(2,:),'r');
    legend('x1','x2');
    grid on;


    figure(6);
    plot([0:h:finaltime-h],E(:),'g-');
    legend('e(t)');

    figure(7);
    plot([0:h:finaltime],XF(1,:),'k-',[0:h:finaltime],XF(2,:),'k.')
    legend('xf(1,:)','xf(2,:)')
    grid on; 
end
figure(8);
plot([0:h:finaltime-h],Bernoulli_Event,'k-')
legend('a(t)')
grid on;

