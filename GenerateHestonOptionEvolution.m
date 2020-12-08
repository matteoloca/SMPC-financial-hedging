function P = GenerateHestonOptionEvolution(strike,w0,N,type,mu,sigma,theta1,k,T,dt,r,y1,barrier,C)
%GENERATEPRICESEVOLUTION Summary of this function goes here
%   Detailed explanation goes here


mdl=gbm(mu,sigma,'StartState',w0);

if (strcmp(type,'barrier')==1)
    OptModel=BarrierOption(T+1,N);
else
    if(strcmp(type,'cliquet')==1)
        
    else
        OptModel=EuropeanOption(T+1,N);
    end
end
mdl.simulate(T,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);

if (strcmp(type,'barrier')==1)
    P=OptModel.OptionPrice(strike,r,barrier);
else
    P=OptModel.OptionPrice(strike,r);
end

end

