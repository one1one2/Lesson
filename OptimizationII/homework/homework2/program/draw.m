plot(log(LBFGS_Hinfo.g)/log(10))
hold on
plot(log(LBFGS_Binfo.g)/log(10))
hold on 
plot(log(LSR1info.g)/log(10))
hold on 
plot(log(PBFGSinfo.g)/log(10))
legend('LBFGS-H', 'LBFGS-B', 'LSR1', 'PBFGS')
xlabel('iter')
ylabel('log |g|')
