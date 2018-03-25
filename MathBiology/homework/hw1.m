E = 0.01;
S = 1;
ES = 0;
P = 0;
E0 = 0.01;

x = [S;ES;P];
tf = 100000;
[t, X] = ode45('differ', [0, tf], x);

S = X(:,1);
ES = X(:,2);
P = X(:,3);
E = E0 - ES;
%plot(t,S,t,P)
%legend('S','P')
%xlabel('t')

P1 = P(2:length(P));
P = P(1:length(P)-1);
t1 = t(2:length(t));
t = t(1:length(P));
res1 = (P1-P)./(t1-t);
res1(length(S)) = 0;
res2 = 0.00001*S./(0.101+S);
plot(S, res1, 'ro', S, res2, 'b')
legend('experiment', 'Michaelis-Menten')
xlabel('[S]');
ylabel('V');
