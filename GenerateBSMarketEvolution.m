function [W,B,P] = GenerateBSMarketEvolution(w0,N,mu,sigma,T,dt,r,opttype,strike,barrier,C,fixDates,allW,allT,sT,overbarrier)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B_tmp=zeros(T,N);
%W=zeros(T+2,N);

if (strcmp(opttype,'barrier')==1)
    OptModel=BarrierOption(T+1,N);
else
    if (strcmp(opttype,'cliquet')==1)
        OptModel=NapoleonCliquet(allT+1,N,fixDates,allW,sT);
    else
        OptModel=EuropeanOption(T+1,N);
    end
end
mdl=gbm(mu,sigma,'StartState',w0);

[tmp,]=mdl.simulate(T,'DeltaTime',dt,'nTrials',N,'Antithetic',false,'Processes',OptModel.Sim); 
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
if P>barrier
    P=0;
end
if overbarrier
    P=0;
end
end

