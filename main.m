clear
[strike,w0,M,Ns,L,T,dt,mu,r,sigma,theta1,...
    k,variance,y1,simtype,opttype,barrier,C]= System_Info;

%P=GenerateHestonOptionEvolution(strike,w0,Ns,type,mu,sigma,theta1,k,T,dt,r,y1,barrier);
%[W,B]=GenerateHestonMarketEvolution(w0,Ns,mu,sigma,theta1,k,T,dt,r,y1);
rng('shuffle');
[P,W,B]=GenerateHestonObjective(strike,w0,Ns,mu,sigma,theta1,k,T,dt,r,y1,opttype,barrier,C);

x=zeros(T,1);
u=zeros(T,1);
e=zeros(T,1); 
p_simulated=zeros(M,T);
w_simulated=zeros(M,T);
x(1)=P(1);
for i = 1:T-1
    [u(i),p_simulated(:,i),~]=LS_Optimization(W(i),P(i),M,Ns,L,x(i),strike,T-i+1,dt,mu,r,sigma,theta1,...
        k,variance,y1,simtype,opttype,barrier,C);
    [x_n,e_n]=next_x(x(i),B(i),u(i),P(i+1),r);
    x(i+1)=x_n;
    e(i+1)=e_n;
end

figure(1)
plot(x)
hold on
plot(P)
p_s= mean(p_simulated)'
plot([P(1); p_s(1:T-1)])
legend('Portfolio wealth','Option price')%,'Option simulated')
title('Portfolio wealth vs option price')

figure(2)
plot(e)
title('Error between portfolio wealth and option price')

mean(e)
mean(abs(e))
max(abs(e))