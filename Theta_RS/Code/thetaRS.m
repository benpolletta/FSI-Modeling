model_eqns = ['dv/dt=I(t)+@current/Cm; Cm=1; v(0)=-65;',...
    '{iNa,iK}; I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=100; toff=225; I_app=10;',...
    'monitor iNa.functions, iK.functions'];
data=SimulateModel(model_eqns,'tspan',[0 1000],'vary',{'aoffset',0:.001:.01; 'boffset',0:.001:.01}); % {'pop1','gd',5:10;'pop1','I_app',10:20});
PlotData(data)

% figure; plot(data.time,data.(data.labels{1}))
% xlabel('time (ms)'); ylabel('membrane potential (mV)'); title('Hodgkin-Huxley neuron w/ K_D channel')