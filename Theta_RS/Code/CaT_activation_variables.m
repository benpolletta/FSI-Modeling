
%%

X = -100:100;

as = 1.6./(1+exp(-.072*(X-65)));
bs = 0.02*(X-51.1)./(exp((X-51.1)/5)-1);
sinf = as./(as+bs);
stau = 1./(as+bs);
% ICa = gCa.*(s.^2).*(X-ECa)

figure

plot(X', sinf')

plot(X', stau')