function [W,B,P] = GenerateHestonMarketEvolution(w0,N,mu,rho1,theta1,k,T,dt,r,y1,omega1,opttype,strike,barrier,C,fixDates,allW,allT,sT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B_tmp=zeros(T,N);
hst=heston(mu,theta1,k,omega1,'StartState',[w0 y1]');
hst.Correlation=[1 rho1; rho1 1];

if (strcmp(opttype,'barrier')==1)
    OptModel=BarrierOption(T+2,N);
else
    if (strcmp(opttype,'cliquet')==1)
        OptModel=NapoleonCliquet(allT+2,N,fixDates,allW,sT);
    else
        OptModel=EuropeanOption(T+2,N);
    end
end
[tmp,]=hst.simulate(T,'DeltaTime',dt,'nTrials',N,'Antithetic',false,'Processes',OptModel.Sim); 
W_tmp=squeeze(tmp(:,1,:));
B_tmp(1:T,:)=W_tmp(2:T+1,:)-((1+r)*W_tmp(1:T,:));
W=mean(W_tmp(2,:));
B=mean(B_tmp(1,:));
B=W-((1+r)*w0);

if (strcmp(opttype,'barrier')==1)
    P=(OptModel.OptionPrice(strike,r,barrier))';
else
    if (strcmp(opttype,'cliquet')==1)
        P=(OptModel.OptionPrice(strike,r,C))';
    else
        P=(OptModel.OptionPrice(strike,r))';
    end
end
P=real(P); %Controllo messo perché ogni tanto risultato complesso
end

