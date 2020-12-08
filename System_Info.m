function [strike,w0,M,Ns,L,T,dt,mu,r,sigma,theta1,k,variance,y1,simtype,opttype,barrier,C] = System_Info
%Genearal data
strike=100;
w0 = 100;
M=50;
Ns=100;
L=1000;
variance=1;
T=24;
dt=1/100;

%Black-scholes model data
mu=.04;
r = .00074102;
sigma = .5;

%Heston model data
k = 1;
ret =.25;
rho1=-0.5;
w1=.3;
theta1=.25;
y1=.25;

%MPC simulation Type: heston or bs
simtype='heston';

%Hedging Type: barrier, european or cliquet
opttype='barrier';

%Barrier data
barrier=120;

%Cliquet data
C = .1;
Tfix=[8,16,24];


end

