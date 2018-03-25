plot(log(LBFGS_Hinfo.f)/log(10))
hold on
plot(log(LBFGS_Binfo.f)/log(10))
hold on 
plot(log(LSR1info.f)/log(10))
hold on 
plot(log(PBFGSinfo.f)/log(10))
legend('LBFGS-H', 'LBFGS-B', 'LSR1', 'PBFGS')
xlabel('iter')
ylabel('log f')

