%% state space
clear; clc;
syms T1 T2
syms m1 m2;
syms a1 a2;
syms I1 I2;
syms L1 L2;
syms c1 s1 c2 s2 c12;
syms theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2;
syms g;
c1=cos(theta1);s1=sin(theta1);c2=cos(theta2);s2=sin(theta2);c12=cos(theta1+theta2)
syms M11 M12 M21 M22 C1 C2 g1 g2
ddtheta=[ddtheta1;ddtheta2]

%dynamic equation
eqn1=[M11*ddtheta1+M12*ddtheta2+C1+g1==T1]
eqn2=[M21*ddtheta1+M22*ddtheta2+C2+g2==T2]
eqn=[eqn1,eqn2]
var=[ddtheta1 ddtheta2];
%solve ddtheta
[ddthe1,ddthe2]=solve(eqn,var)
pretty(simplify(ddthe1))             %ddtheta1
pretty(simplify(ddthe2))           %ddtheta2
%% linearization
% equilibrium
clear; clc;
L1=1;L2=1;a1=0.5;a2=0.5;I1=0.25;I2=0.25;m1=1;m2=1;g=9.8;
T1=10;T2=10
% syms T1 T2
% T1=10;T2=10
% syms m1 m2;
% syms a1 a2;
% syms I1 I2;
% syms L1 L2;
% syms c1 s1 c2 s2 c12;
syms theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2;
% syms g;
% theta1=deg2rad(30)
% theta2=deg2rad(60)
% dtheta1=0
% dtheta2=0
c1=cos(theta1);s1=sin(theta1);c2=cos(theta2);s2=sin(theta2);c12=cos(theta1+theta2)
M11=m1*L1^2+I1+m2*(a1^2+L2^2+2*a1*a2*c2)+I2;
M12=m2*(L2^2+a1*L2*c2)+I2;
M21=m2*(L2^2+a1*L2*c2)+I2;
M22=m2*L2^2+I2;
C1=-m2*a1*L2*s2*(2*dtheta1*dtheta2+dtheta2^2);
C2=m2*a1*L2*s2*dtheta1^2;
g1=m1*g*L1*c1+m2*g*(a1*c1+L2*c12);
g2=m2*g*L2*c12;
% To linearized the system, find equalibrium point, where dtheat1 and dtheta2=0
% state space=0, using taylor expansion
equal1=[taylor(M12*T2 - M22*T1 - M12*g2 + M22*g1,[theta1,theta2,dtheta1,dtheta2],'order',4)==0]
equal2=[taylor(M11*T2 - M21*T1 - M11*g2 + M21*g1,[theta1,theta2,dtheta1,dtheta2],'order',4)==0]
equal=[equal1,equal2];
var=[theta1 theta2];
% solve equlibrium point:theta1 theta2
[the1,the2]=solve(equal,var) 

%Calculate A ,B matrix by taking Jacobian of state space,and plug in equalibrium point calculated above
syms T1 T2
ddthe1=-(M12*T2 - M22*T1 - M12*g2 + M22*g1 - C2*M12 + C1*M22)/(M11*M22 - M12*M21)
ddthe2=(M11*T2 - M21*T1 - M11*g2 + M21*g1 - C2*M11 + C1*M21)/(M11*M22 - M12*M21)
f1=dtheta1;        
f2=dtheta2;
f3=ddthe1;
f4=ddthe2;
A=jacobian([f1;f2;f3;f4],[theta1 theta2 dtheta1 dtheta2])
B=jacobian([f1;f2;f3;f4],[T1 T2])
%plug in equilibrium point
theta1_eqn=deg2rad(the1(4))
theta2_eqn=deg2rad(the2(4))
dtheta1_eqn=0
dtheta2_eqn=0
A_lin=subs(A,{theta1,theta2,dtheta1,dtheta2,T1,T2},{theta1_eqn,theta2_eqn,dtheta1_eqn,dtheta2_eqn,10,10}) 
B_lin=subs(B,{theta1,theta2,dtheta1,dtheta2,T1,T2},{theta1_eqn,theta2_eqn,dtheta1_eqn,dtheta2_eqn,10,10})
C_lin=[1 0 0 0;
   0 1 0 0];
D_lin=[];
sys=ss(double(vpa(A_lin,3)),double(vpa(B_lin,3)),C_lin,D_lin)

%discrete model
Ts=0.05
sysD=c2d(sys,Ts)        %continuous to discrete state space
tf=1
td = Ts:Ts:tf;
x = zeros(4,length(td));
y = zeros(2,length(td));
u1=10;u2=10
u=[u1;u2]
for k = 2:length(x)
x(:,k) = sysD.A * x(:,k-1) + sysD.B*u;
y(:,k) = sysD.C*x(:,k);
end
figure(3)
stem(td, y(1,:),'r')
xlabel('time(s)');
ylabel('rad')
title('theta1');
figure(4)
stem(td, y(2,:),'b')
xlabel('time(s)');
ylabel('rad')
title('theta2');

