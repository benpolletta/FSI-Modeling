function data = Gutfreund_fake_dendrite(Iconst, I_app, varargin)

vary_label = ''; vary_cell = cell(length(varargin)/2, 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        vary_label = [vary_label, sprintf('_%s_%fto%f', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        vary_cell(a, :) = {'pop1', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

model_eqns = ['dv/dt=Iconst+I(t)+@current/Cm; Cm=.25; v(0)=-65;',...
    sprintf('{iNaP,iKs,iKDRG,iNaG,gleak,faux_Gutfreund_dendrite}; halfKs=-35; halfNaP=-50; gNaP=0.025; Iconst=%f;', Iconst),...
    'offset=0; Koffset=offset; Noffset=offset; gdenom=1; gKDR=5/gdenom; gNa=12.5/gdenom; gl=0;',... % gdenom=3; gKDR=5/gdenom; gNa=12.5/gdenom; gl=0;'
    'tau_fast=7; tauh=tau_fast; taum=tau_fast;',...
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell);
    
else
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000]); %{'pop1','gKs',.04:.002:.06 % ;'pop1','gNaP',.010:.001:.020}); % {'pop1','gd',5:10;'pop1','I_app',10:20});

end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('gutfreund_fd_Iconst%f_Iapp%f_', Iconst, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end