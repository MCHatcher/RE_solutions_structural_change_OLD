% Inflation target plotter (credibility example)

for t=1:T_sim-1

        pi(t) = X_stack(1,t); y(t) = X_stack(2,t); int(t) = X_stack(3,t); 
        
end

%set(0,'DefaultLineLineWidth',1.5)
subplot(2,3,1), hold on, plot(Periods, 4*100*pi,'--k'), title('Inflation'), xlabel('Periods'), ylabel('% (annualized)'), axis([-inf,inf,-inf,inf]), %axis([1,16,2,6])
subplot(2,3,2), hold on, plot(Periods, 100*y,'--k'), title('Output'), xlabel('Periods'), ylabel('%'), axis([-inf,inf,-inf,inf]), %axis([1,16,-0.2,0.6])
subplot(2,3,3), hold on, plot(Periods, 4*100*int,'--k'), title('Nominal interest rate'), axis([-inf,inf,-inf,inf])  %axis([1,16,5,9.01]),
xlabel('Periods'), ylabel('% (annualized)')