function data = Gutfreund(I_const, I_app, varargin)

vary_label = ''; vary_cell = {};

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        vary_label = [vary_label, sprintf('_%s_%fto%f', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        vary_cell = {vary_cell{:}; 'pop1', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

model_eqns = ['dv/dt=I_const+I(t)+@current/Cm; Cm=.25; v(0)=-65;',...
    sprintf('{iNaP,iKs}; halfKs = -60; I_const=%f;', I_const),...
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=1500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell);
    
else
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000]); %{'pop1','gKs',.04:.002:.06 % ;'pop1','gNaP',.010:.001:.020}); % {'pop1','gd',5:10;'pop1','I_app',10:20});

end
    
PlotData(data)

save_as_pdf(gcf, [sprintf('gutfreund_Iconst%f_Iapp%f_', I_const, I_app), vary_label])