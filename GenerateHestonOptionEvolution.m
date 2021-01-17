function [P,W] = GenerateHestonOptionEvolution(strike,w0,N,type,mu,rho1,theta1,k,T,dt,r,y1,omega1,barrier,C,fixDates)
%GENERATEPRICESEVOLUTION Summary of this function goes here
%   Detailed explanation goes here

hst=heston(mu,theta1,k,omega1,'StartState',[w0 y1]');
hst.Correlation=[1 rho1; rho1 1];

if (strcmp(type,'barrier')==1)
    OptModel=BarrierOption(T+2,N);
else
    if(strcmp(type,'cliquet')==1)
        OptModel=NapoleonCliquet(T+2,N,fixDates);
    else
        OptModel=EuropeanOption(T+2,N);
    end
end
hst.simulate(T+1,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);

if (strcmp(type,'barrier')==1)
    P=OptModel.OptionPrice(strike,r,barrier);
else
    if (strcmp(type,'cliquet')==1)
        P=(OptModel.OptionPrice(strike,r,C))';
    else
        P=OptModel.OptionPrice(strike,r);
    end
end
W_t=OptModel.Prices();
W=W_t(:,2);
end

