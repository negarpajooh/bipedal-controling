clear;
clc;
close all
format long
global params State setpoint_choose Teta As Bs Cs Ds
%% CTE
m1=0.5;m2=1;m3=1;m4=1;m6=1;
L1d=1;L2d=1;L3d=1;L4d=1;L6d=1;
L1=1;L2=1;L4=1;L5=1;L6=1;
I1=1;I2=1;I3=1;I4=1;I6=1;
L1u=1;L2u=1;L4u=1;L5u=1;L6u=1;
b1=1;b2=1;b3=1;mp2=1;g=10;
%%
params=[m1,m2,m3,m4,m6,L1d,L2d,L3d,L4d,L6d,L1,L2,L4,L5,L6,I1,I2,I3,I4,I6,L1u,L2u,L4u,L5u,L6u,b1,b2,b3,mp2,g];
%%
D=zeros(5);D(1,1)=1;D(1,2)=1;D(2,2)=-1;D(2,3)=-1;D(3,4)=1;D(4,4)=-1;D(4,5)=-1;D(5,3)=1;D(5,5)=1;
%% Initial Condotion
State=questdlg('Plz Selection State Bipedal?','Initial Value','1','2',1)
switch State
    case '1'
        Teta=[0,0,0,0,0];
        setpoint_choose=3;
    case '2'
        Teta=[pi/6,pi/4,0,0,0];
        setpoint_choose=0;
end
TetaDot=[0,0,0,0,0];
%% M G C
M=Mass(params,Teta,TetaDot);
G=stiffness(params);
C=Damping(params,Teta,TetaDot);
[mi,~]=size(M);
%% State Space
As=[zeros(mi),eye(mi);-G*pinv(M),-C*pinv(M)];
[na,~]=size(As);
Bs=[zeros(.5*na);D*pinv(M)];
Cs=eye(na);
Ds=zeros(10,5);
%%
Sys=ss(As,Bs,Cs,Ds);
save('Sys')
mpctool
