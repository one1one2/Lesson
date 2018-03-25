V = [0.15, 0.34, 0.35, 0.42, 0.48, 0.60, ...
     0.52, 0.63, 0.63, 0.63, 0.60, 0.66, ...
     0.69, 0.63, 0.73, 0.69, 0.74, 0.77, ...
     0.72, 0.75 ];
S = [0.10, 0.15, 0.19, 0.24, 0.29, 0.34, ...
     0.38, 0.43, 0.48, 0.53, 0.57, 0.62, ...
     0.67, 0.72, 0.76, 0.81, 0.86, 0.91, ...
     0.95, 1.00];

g = fittype('Vmax*x/(Km+x)');
res1 = fit(S', V', g)
Vmax = res1.Vmax
Km = res1.Km 
res2 = fit(1./S', 1./V', 'poly1');
Vmax = 1/res2.p2
Km = p1/p2

plot(S, V, 'ko', S, res1(S), 'r', S, Vmax*S./(Km+S), 'b')     
legend('initial', 'Mochaelis', 'Lineweaver')
xlabel('[S]')
ylabel('V')


