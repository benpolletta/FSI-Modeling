g_struct.gK = 225*[2; 0]; g_struct.gNa = 112.5*[2; 0]; g_struct.gL = 0.25*[2; 1/10];

[Vs,Vd,~,~,~,~,~,t] = Golomb_2007_RK45_w_dendrite(1,27,2000,g_struct,0);
figure; plot(t,Vs), hold on, plot(t,Vd,'k')
[Vs,Vd,~,~,~,~,~,t] = Golomb_2007_RK45_w_dendrite(1,27,2000,struct(),0);
hold on, plot(t,Vs,'c'), hold on, plot(t,Vd,'m')
[Vs,Vd,~,~,~,~,~,t] = Golomb_2007_RK45_w_dendrite(1,27,2000,struct('g_sd',0),0);
hold on, plot(t,Vs,'g'), hold on, plot(t,Vd,'r')
title('Frequency: No Dendrite > Passive Dendrite > Active Dendrite')

save_as_pdf(gcf,'passive_v_active_v_none')