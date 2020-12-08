function [P,W,B] = GenerateHestonObjective(strike,w0,N,mu,sigma,theta1,k,T,dt,r,y1,type,barrier,C)
%GENERATEHESTONOBJECTIVE Summary of this function goes here
%   Detailed explanation goes here
B=zeros(T-1,1);
P=zeros(T,1);
W=zeros(T,1);
%First run
if (strcmp(type,'barrier')==1)
    OptModel=BarrierOption(T+1,N);
else
    OptModel=EuropeanOption(T+1,N);
end
W(1)=w0;

hst=heston(mu,k,theta1,sigma,'StartState',[w0 y1]');

hst.simulate(T,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);
if (strcmp(type,'barrier')==1)
    P(1)=(OptModel.OptionPrice(strike,r,barrier))';
else
    P(1)=(OptModel.OptionPrice(strike,r))';
end
% W_tmp=(OptModel.Prices());
% W=mean(W_tmp);
% 
% if (strcmp(type,'barrier')==1)
%     P=(OptModel.AllOptionPrices(strike,r,barrier,dt));
% else
%     P=(OptModel.AllOptionPrices(strike,r,dt));
% end

% %All others
for i=1:T-1
    if (strcmp(type,'barrier')==1)
        OptModel=BarrierOption(T+1-i,N);
    else
        OptModel=EuropeanOption(T+1-i,N);
    end
    hst=heston(mu,k,theta1,sigma,'StartState',[W(i) y1]');
    hst.simulate(T-i,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);
    if (strcmp(type,'barrier')==1)
        P(i+1)=(OptModel.OptionPrice(strike,r,barrier))';
    else
        P(i+1)=(OptModel.OptionPrice(strike,r))';
    end
    W_tmp=(OptModel.Prices());
    W(i+1)=mean(W_tmp(:,2));
end
%P(T+1)=max(W(T)-strike,0);

B(1:T-1)=W(2:T)-((1+r)*W(1:T-1));

end








