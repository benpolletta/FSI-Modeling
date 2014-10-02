input = rand(10, ceil(100/.005));

input = input < 1000*.005/1000;

for i = 1:10, input(i, :) = conv(double(input(i,:)), ones(1,200), 'same'); end

[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(10, input*70, 100, [], zeros(10), .3*(rand(10) > .4));

plot(Vs')