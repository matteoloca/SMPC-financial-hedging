function [W,B] = GenerateHestonMarketEvolution(w0,N,mu,sigma,theta1,k,T,dt,r,y1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B_tmp=zeros(T,N);
%W=zeros(T+2,N);

hst=heston(mu,k,theta1,sigma,'StartState',[w0 y1]');

[tmp,time]=hst.simulate(T,'DeltaTime',dt,'nTrials',N,'Antithetic',false); 
W_tmp=squeeze(tmp(:,1,:));
B_tmp(1:T,:)=W_tmp(2:T+1,:)-((1+r)*W_tmp(1:T,:));
W=mean(W_tmp(2,:));
B=mean(B_tmp(1,:));
B=W-((1+r)*w0);
end

