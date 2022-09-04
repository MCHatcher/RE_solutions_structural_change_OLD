% Inflation target plotter

for t=1:T_sim

        pi(t) = X_stack(1,t); y(t) = X_stack(2,t); int(t) = X_stack(3,t);
        
        pi_orig(t) = X_stack_orig(1,t); y_orig(t) = X_stack_orig(2,t); 
        int_orig(t) = X_stack_orig(3,t);
        
        pi_fin(t) = X_stack_fin(1,t); y_fin(t) = X_stack_fin(2,t);
        int_fin(t) = X_stack_fin(3,t);
        
end 

set(0,'DefaultLineLineWidth',2)
hold on, 
subplot(1,3,1), hold on, plot(Periods, 4*100*pi, 'k'), hold on, plot(Periods, 4*100*pi_orig,'--k'), 
hold on, plot(Periods, 4*100*pi_fin,'.-.k'), title('Inflation'), xlabel('Periods'), ylabel('% (annualized)'), axis([1,inf,2,6]) %axis([-inf,inf,-inf,inf])
subplot(1,3,2), hold on, plot(Periods, 100*y,'k'), hold on, plot(Periods, 100*y_orig, '--k'), title('Output'),
xlabel('Periods'), ylabel('%'), axis([1,inf,-0.2,0.6]), %axis([-inf,inf,-inf,inf]) 
subplot(1,3,3), hold on, plot(Periods, 4*100*int,'k'), hold on, plot(Periods, 4*100*int_orig,'--k'), 
hold on, plot(Periods, 4*100*int_fin,'.-.k'), title('Nominal interest rate'),  axis([1,inf,5,9.01]),  %axis([-inf,inf,-inf,inf])
xlabel('Periods'), ylabel('% (annualized)'), hold on