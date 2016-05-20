function data = Gutfreund(I_const, I_app, varargin)

vary_label = ''; vary_cell = cell(length(varargin)/2, 3);

if ~isempty(varargin)
    
    for a = 1:(length(varargin)/2)
        
        vary_label = [vary_label, sprintf('_%s_%fto%f', varargin{2*a - 1}, varargin{2*a}(1), varargin{2*a}(end))];
        
        vary_cell(a, :) = {'pop1', varargin{2*a - 1}, varargin{2*a}};
        
    end
    
end

model_eqns = ['dv/dt=I_const+I(t)+@current/Cm; Cm=.25; v(0)=-65;',...
    sprintf('{iNaP,iKs,iKDRG,iNaG,gleak}; halfKs=-60; halfNaP=-60; gNaP=0.0125; I_const=%f;', I_const),...
    'offset=-16; Koffset=offset; Noffset=offset; gKDR=5/3; gNaG=12.5/3; gl=.025/2;',...
    sprintf('I(t)=I_app*(ton<t&t<toff)*(1+rand*.25); ton=500; toff=3500; I_app=%f;', I_app),...
    'monitor functions'];

if ~isempty(varargin)
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000], 'vary', vary_cell);
    
else
    
    data=SimulateModel(model_eqns, 'tspan', [0 4000]); %{'pop1','gKs',.04:.002:.06 % ;'pop1','gNaP',.010:.001:.020}); % {'pop1','gd',5:10;'pop1','I_app',10:20});

end
    
try 
    
    PlotData(data)

    save_as_pdf(gcf, [sprintf('gutfreund_Iconst%f_Iapp%f_', I_const, I_app), vary_label])

catch error
    
    display('PlotData failed:')
    
    display(error)

end