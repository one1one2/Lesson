function y = differ(t, x)
  k1p = 1;
  k1m = 0.1;
  k2 = 0.001;
  E0 = 0.01;
  S = x(1);
  ES = x(2);
  P = x(3);
  E = E0 - ES;
  dS = -k1p*E*S + k1m*ES;
  dES = k1p*E*S - (k1m+k2)*ES;
  dP = k2*ES;
  y = [dS;dES;dP];
end 

