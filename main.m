clear
[strike,w0,M,Ns,L,T,dt,mu,r,sigma,theta1,...
    k,y1,omega1,rho1,simtype,opttype,barrier,C,fixDates]= System_Info;

rng('shuffle');
%rng('default')
if(strcmp(opttype,'european')==1)
    N=1000;
else
    N=1000;
end
[P,W,B,variance]=GenerateHestonObjective(strike,w0,N,mu,theta1,k,T,dt,r,y1,rho1,omega1,opttype,barrier,C,fixDates);

x=zeros(T+1,1);
if (strcmp(opttype,'barrier')==1 || strcmp(opttype,'cliquet')==1)
    u=zeros(T,2);
else
    u=zeros(T,1);
end
e=zeros(T+1,1);
p_simulated=zeros(M,T+1);
w_simulated=zeros(M,T+1);
err_funct=zeros(T,1);
e_medio=zeros(T,1);
x(1)=P(1);
for i = 1:T
    [u(i,:),p_simulated(:,i),~,err_funct(i),e_medio(i)]=LS_Optimization(W(i),w0,P(i),M,Ns,L,x(i),strike,T+2-i,i,dt,mu,r,sigma,theta1,...
        k,variance(i),omega1,rho1,simtype,opttype,barrier,C,fixDates,W(1:i),T);
    [x_n,e_n]=next_x(x(i),B(i,:)',u(i,:)',P(i+1),r);
    x(i+1)=x_n;
    e(i+1)=e_n;
end

figure(1)
plot([1:25],x);%,'color','g')
hold on
plot([1:25],P);%,'color','r');
p_s= mean(p_simulated)';
%plot([P(1); p_s(1:T)])
legend('Portfolio wealth','Option price')%,'Option simulated')
title('Portfolio wealth vs option price')
ylabel('Portfolio wealth')%, option price')
xlabel('Time(Weeks)')

figure(2)
%ef_m=mean(err_funct);
plot(e);
hold on
%plot(err_funct)
%plot([0; e_medio])
title('Error between portfolio wealth and option price')
xlabel('Time(Weeks)')
ylabel('Error')
% figure(3)
% plot(u)
% hold on
% plot(W)
figure(3)
plot(p_simulated')

mean(e)
mean(abs(e))
max(abs(e))