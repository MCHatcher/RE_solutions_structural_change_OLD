% Forward guidance plotter

for t=1:T_sim

        pi(t) = X_stack(1,t); y(t) = X_stack(2,t); int(t) = X_stack(3,t);
        
        pi_fin(t) = X_stack_fin(1,t); y_fin(t) = X_stack_fin(2,t);
        int_fin(t) = X_stack_fin(3,t);
        
end 

set(0,'DefaultLineLineWidth',1.5)
hold on, 
subplot(1,3,1), hold on, plot(Periods, 4*100*pi, '--'), hold on, %plot(Periods, 4*100*pi_fin,'--k'), 
title('Inflation'), xlabel('Periods'), ylabel('% (annualized)'), axis([-inf,inf,-inf,inf]) %axis([-inf,inf,-inf,inf])
subplot(1,3,2), hold on, plot(Periods, 100*y,'--'), hold on, %plot(Periods, 100*y_fin, '--k'), 
title('Output'),xlabel('Periods'), ylabel('%'), axis([-inf,inf,-inf,inf]), %axis([-inf,inf,-inf,inf]) 
subplot(1,3,3), hold on, plot(Periods, 4*100*int,'--'), hold on, %plot(Periods, 4*100*int_fin,'--k'), 
title('Nominal interest rate'),  axis([-inf,inf,-inf,inf]),  %axis([-inf,inf,-inf,inf])
xlabel('Periods'), ylabel('% (annualized)'), hold on