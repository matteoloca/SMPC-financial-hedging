function [u_n,p_sim,w_sim] = LS_Optimization(w0,p0,M,Ns,L,x,strike,T,dt,mu,r,sigma,theta1,...
        k,variance,y1,simtype,opttype,barrier,C)
%LS_OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here
if (strcmp(opttype,'barrier')==1)
    w_sim=zeros(2,M);
    b_sim=zeros(2,M);
else
    w_sim=zeros(1,M);
    b_sim=zeros(1,M);
end

p_sim=zeros(1,M);

for i=1 : M
    if (strcmp(simtype,'heston')==1)
        [w_sim(1,i),b_sim(i)]=GenerateHestonMarketEvolution(w0,Ns,mu,sigma,theta1,k,T,dt,r,y1);
        p_sim(i)=real(GenerateHestonOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));
    else
        [w_sim(1,i),b_sim(i)]=GenerateBSMarketEvolution(w0,Ns,mu,sigma,theta1,k,T,dt,r,y1);
        p_sim(i)=real(GenerateBSOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));

    end
    if (strcmp(opttype,'barrier')==1)
        if (strcmp(simtype,'heston')==1)
            w_sim(2,i)=real(GenerateHestonOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));
            b_sim(2,i)=w_sim(2,i)-((1+r)*p0);
        else
            w_sim(2,i)=real(GenerateBSOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));
            b_sim(2,i)=w_sim(2,i)-((1+r)*p0);
        end
    end
    
    end

%w_sim(i)=w_sim(2,:);
%b_sim(i)=b_sim(1,:);

%solve the LS problem
u=optimvar('u');
prob=optimproblem;
[x_n,e]=next_x(x,b_sim,u,p_sim,r);
%f=sum((1/M)*((e-(sum(e)/M)).^2));
f=sum((1/M)*((e-mean(e)).^2));
%f=var((e-mean(e)).^2);
%f=mean(e);
expr=f;
prob.Objective=expr;
options = optimoptions('lsqlin','Display','off');
sol=solve(prob,'Options',options);
u_n=sol.u;

[x_n,e]=next_x(x,b_sim,u_n,p_sim,r);
expr=sum((1/M).*((e-(sum(e)/M)).^2));

end

