function [strike,w0,M,Ns,L,T,dt,mu,r,sigma,theta1,k,y1,omega1,rho1,simtype,opttype,barrier,C,Tfix] = System_Info
%Genearal data
strike=100;
w0 = 100;
M=100;
Ns=50;
L=50;
variance=1;
T=24;
dt=1/52;

%Black-scholes model data
mu=.04;
r = .00074102;
sigma = .5;

%Heston model data
k = 1;
rho1=-0.5;
omega1=.3;
theta1=.25;
y1=.25;

%MPC simulation Type: heston or bs
simtype='heston';

%Hedging Type: barrier, european or cliquet
opttype='cliquet';

%Barrier data
barrier=120;

%Cliquet data
C = .1;
Tfix=[1,9,17,25];


end

