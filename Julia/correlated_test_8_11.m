firing_rate = zeros(10,10);
firing_avg = zeros(10,1);
gap_conductance = 0.3;
rhythmicity = zeros(10,10,20000,10);
for k = 1:10
	for j = 1:10
		input = rand(10, ceil(100/.005));
		shared_input = rand(1, ceil(100/.005));
		input = input < 1000*.005/1000;  
		shared_input = shared_input < 1000*.005/1000;
		for i = 1:10, input(i, :) = conv(double(input(i,:)), ones(1,200), 'same'); end
		shared_input(1, :) = conv(double(input(1,:)), ones(1,200), 'same');
		fraction_shared = 0.1*k;
		for i = 1:10, input(i, :) = fraction_shared*shared_input + (1-fraction_shared)*input(i,:); end
		[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(10, input*70, 100, [], zeros(10), gap_conductance*(rand(10) > .4));
		rhythmicity(k,j,:,:) = abs(fft(Vs'));
		firing_rate(k,j) = 0;
		for a = 1:10
			for b = 1:19999
				if Vs(a,b) < 0 && Vs(a,b+1) > 0
        			firing_rate(k,j)= firing_rate(k,j) + 1;
				end
			end
        end
        firing_rate(k,j) = firing_rate(k,j)/10
    end
    firing_avg = sum(firing_rate,2)/10
end

for i = 1:10
    for j = 1:10
        subplot(10,1,i)
        plot(squeeze(rhythmicity(i,j,:,:)))
    end
end