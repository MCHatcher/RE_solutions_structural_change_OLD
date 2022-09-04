%Pension reform plotter

for t=1:T_sim
    
        cy(t) = X_stack(1,t); r(t) = X_stack(2,t); w(t) = X_stack(3,t);
        k(t) = X_stack(4,t); co(t) = X_stack(5,t);
        klev(t) = exp(k(t) + log(kstar));
        
        cyu(t) = Xu_stack(1,t); ru(t) = Xu_stack(2,t); wu(t) = Xu_stack(3,t);
        ku(t) = Xu_stack(4,t); cou(t) = Xu_stack(5,t);
       
end 

set(0,'DefaultLineLineWidth',1.5)
figure(1)
subplot(1,3,1), plot(Periods(1:12),cy(1:12),'k'), hold on, plot(Periods(1:12),cyu(1:12),'--k'), title('Consumption (young): $\hat{c}_{t,y}$','interpreter','latex'), xlabel('Time,t')
subplot(1,3,2), plot(Periods(1:12),co(1:12),'k'),hold on, plot(Periods(1:12),cou(1:12),'--k'), title('Consumption (old): $\hat{c}_{t,o}$', 'interpreter','latex'), xlabel('Time,t')
subplot(1,3,3), hold on, plot(Periods(1:12), k(1:12),'k'), hold on, plot(Periods(1:12),ku(1:12),'--k'), title('Capital: $\hat{k}_{t+1}$','interpreter','latex'), xlabel('Time,t')
