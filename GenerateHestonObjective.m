function [P,W,B,var] = GenerateHestonObjective(strike,w0,N,mu,theta1,k,T,dt,r,y1,omega1,rho1,type,barrier,C,fixDates)

B=zeros(T,1);
P=zeros(T+1,1);
W=zeros(T+1,1);
var=zeros(T+1,1);
%First run
if (strcmp(type,'barrier')==1)
    OptModel=BarrierOption(T+2,N);
else
    if (strcmp(type,'cliquet')==1)
        OptModel=NapoleonCliquet(T+2,N,fixDates,0,1);
    else
        OptModel=EuropeanOption(T+2,N);
    end
end
W(1)=w0;
var(1)=y1;

hst=heston(mu,theta1,k,omega1,'StartState',[w0 y1]');
hst.Correlation=[1 rho1; rho1 1];

hst.simulate(T+1,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);
if (strcmp(type,'barrier')==1)
    P(1)=(OptModel.OptionPrice(strike,r,barrier))';

else
    if (strcmp(type,'cliquet')==1)
        P(1)=(OptModel.OptionPrice(strike,r,C))';
    else
        P(1)=(OptModel.OptionPrice(strike,r))';
    end
end
W_tmp=(OptModel.Prices());
W(2)=mean(W_tmp(:,2));
Var_tmp=(OptModel.Variance());
var(2)=mean(Var_tmp(:,2));

% W_tmp=(OptModel.Prices());
% W=mean(W_tmp);
% 
% if (strcmp(type,'barrier')==1)
%     P=(OptModel.AllOptionPrices(strike,r,barrier,dt));
% else
%     P=(OptModel.AllOptionPrices(strike,r,dt));
% end

% %All others
for i=2:T
    if (strcmp(type,'barrier')==1)
        OptModel=BarrierOption(T+2-i,N);
    else
        if (strcmp(type,'cliquet')==1)
            OptModel=NapoleonCliquet(T+2,N,fixDates,W(1:i),i);
        else
            OptModel=EuropeanOption(T+2-i,N);
        end
    end
    hst=heston(mu,theta1,k,omega1,'StartState',[W(i) var(i)]');
    hst.Correlation=[1 rho1; rho1 1];
    hst.simByTransition(T+2-i,'DeltaTime',dt,'nTrials',N,'Processes',OptModel.Sim);
    if (strcmp(type,'barrier')==1)
        P(i)=(OptModel.OptionPrice(strike,r,barrier))';
    else
        if (strcmp(type,'cliquet')==1)
            P(i)=(OptModel.OptionPrice(strike,r,C))';
        else
            P(i)=(OptModel.OptionPrice(strike,r))';
        end
    end
    W_tmp=(OptModel.Prices());
    W(i+1)=mean(W_tmp(:,2));
    Var_tmp=(OptModel.Variance());
    var(i+1)=mean(Var_tmp(:,2));
    
    if (strcmp(type,'barrier')==1)
        if W(i+1)>barrier
            P(i)=0;
        end
    end
end

if (strcmp(type,'cliquet')==1)
    times = size(fixDates);
    times=times(2);
    %po=zeros(T+1);

    po=((W(fixDates(2:times))-W(fixDates(1:times-1)))./(W(fixDates(1:times-1))));
    value = exp(-r*T)*max(C+min(po),0);
    P(T+1)=value;
    
else
    P(T+1)=max(W(T+1)-strike,0);
    if (strcmp(type,'barrier')==1)
        if W(T+1)>barrier
            P(T+1)=0;
        end
    end
end
P=real(P);

if (strcmp(type,'barrier')==1 || strcmp(type,'cliquet')==1)
    B=zeros(T-1,2);
    B(1:T,2)=P(2:T+1)-((1+r)*P(1:T));
else
    B=zeros(T-1,1);
end
B(1:T,1)=W(2:T+1)-((1+r)*W(1:T));


end








