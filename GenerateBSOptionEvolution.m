function P = GenerateBSOptionEvolution(strike,w0,N,type,mu,sigma,T,dt,r,barrier,C,fixDates)
%GENERATEPRICESEVOLUTION Summary of this function goes here
%   Detailed explanation goes here


mdl=gbm(mu,sigma,'StartState',w0);

if (strcmp(type,'barrier')==1)
    OptModel=BarrierOption(T+1,N);
else
    if(strcmp(type,'cliquet')==1)
        OptModel=NapoleonCliquet(T+1,N,fixDates);
    else
        OptModel=EuropeanOption(T+1,N);
    end
end
mdl.simulate(T,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);

if (strcmp(type,'barrier')==1)
    P=OptModel.OptionPrice(strike,r,barrier);
else
    if (strcmp(type,'cliquet')==1)
        P=(OptModel.OptionPrice(strike,r,C))';
    else
        P=OptModel.OptionPrice(strike,r);
    end
end

end

