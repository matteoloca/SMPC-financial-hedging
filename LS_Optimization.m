function [u_n,p_sim,w_sim,val_expr,e_m] = LS_Optimization(w_act,w0,p0,M,Ns,L,x,strike,T,act_t,dt,mu,r,sigma,theta1,...
        k,variance,omega1,rho1,simtype,opttype,barrier,C,fixDates,allW,allT)

if (strcmp(opttype,'barrier')==1 || strcmp(opttype,'cliquet')==1)
    w_sim=zeros(2,M);
    b_sim=zeros(2,M);
    u=optimvar('u',2,1);
else
    w_sim=zeros(1,M);
    b_sim=zeros(1,M);
    u=optimvar('u',1,1);
end

p_sim=zeros(1,M);

overbarrier=w_act>barrier;
for i=1 : M
    
    if (strcmp(simtype,'heston')==1)
        [w_sim(1,i),b_sim(1,i),p_sim(i)]=GenerateHestonMarketEvolution(w_act,Ns,mu,rho1,theta1,k,T,dt,r,variance,omega1,opttype,strike,barrier,C,fixDates,allW,allT,act_t,overbarrier);
        %p_sim(i)=real(GenerateHestonOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));
    else
        [w_sim(1,i),b_sim(1,i),p_sim(i)]=GenerateBSMarketEvolution(w_act,Ns,mu,sigma,T,dt,r,opttype,strike,barrier,C,fixDates,allW,allT,act_t,overbarrier);
        %p_sim(i)=real(GenerateBSOptionEvolution(strike,w0,Ns,opttype,mu,sigma,theta1,k,T,dt,r,y1,barrier));

    end
    if (strcmp(opttype,'barrier')==1 || strcmp(opttype,'cliquet')==1)
        strike_t=w0*((1+r)^(T-act_t));
        if (strcmp(simtype,'heston')==1)
            w_sim(2,i)=real(GenerateHestonOptionEvolution(strike_t,w_act,L,'European',mu,rho1,theta1,k,1,dt,r,variance,omega1,barrier,C,fixDates));
            b_sim(2,i)=w_sim(2,i)-((1+r)*p0);
        else
            w_sim(2,i)=real(GenerateBSOptionEvolution(strike_t,w_act,L,'European',mu,sigma,1,dt,r,barrier,C,fixDates));
            b_sim(2,i)=w_sim(2,i)-((1+r)*p0);
        end
    end
end

%w_sim(i)=w_sim(2,:);
%b_sim(i)=b_sim(1,:);

%solve the LS problem

prob=optimproblem('ObjectiveSense','min');
[x_n,e_n]=next_x(x,b_sim,u,p_sim,r);
%f=sum((1/M)*((e_n-(sum(e_n)/M)).^2));
%f=sum((1/M)*((e-mean(e)).^2));
%f=var((e-mean(e)).^2);
%f=mean(e);
%f=var(e_n);
%f=e;
%if (strcmp(opttype,'barrier')==1 || strcmp(opttype,'cliquet')==1)
%     function f=vectorObjective(u)
%     %function f=scalarobjective(u)
     f=0;
 for k=1:M
     tmp=(e_n(k)-mean(e_n))^2;
     %[x_n,e_n]=next_x(x,b_sim,u,p_sim,r);
     %f=f+sum((1/M)*((e_n-(sum(e_n)/M)).^2));
     %end
%end
f=f+ (tmp/M);
end
% f=sum((e_n - sum(e_n,1)./M).^2, 1) ./ (M-1);
 expr=f;
 prob.Objective=expr;
 options = optimoptions('QUADPROG','Display','off');
 sol=solve(prob,'Options',options);
 u_n=sol.u;
%options = optimset('Display','iter');
%u_n=fminbnd(@vectorObjective,-100,100,options)
if (size(u_n)==0)
    u_n=0;
end

[x_n,e]=next_x(x,b_sim,u_n,p_sim,r);
val_expr=0;
for k=1:M
    tmp=(e(k)-mean(e))^2;
    val_expr=val_expr+ (tmp/M);
end
e_m=mean(e);
end

