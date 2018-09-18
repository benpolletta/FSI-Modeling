function results = reference_comparison(data, varargin)

center_flag = 1; norm_flag = 2; distance_flag = 2;

if ~isempty(varargin)
    
    for arg = 1:(length(varargin)/2)
        
        if strcmp(varargin{2*arg - 1}, 'variable')
            
            variable = varargin{2*arg};
            
        elseif strcmp(varargin{2*arg - 1}, 'ref_data')
            
            ref_data = varargin{2*arg};
            
        elseif strcmp(varargin{2*arg - 1}, 'center_flag')
            
            center_flag = varargin{2*arg};
            
        elseif strcmp(varargin{2*arg - 1}, 'norm_flag')
            
            norm_flag = varargin{2*arg};
            
        elseif strcmp(varargin{2*arg - 1}, 'distance_flag')
            
            norm_flag = varargin{2*arg};
            
        end
        
    end
    
end

if ~exist('variable', 'var'), variable = 'pop1_v'; end

v = getfield(data, variable);

if size(v, 1) ~= length(v), v = v'; end

if ~exist('ref_data', 'var'), ref_data = zeros(size(v)); end

if size(ref_data, 1) ~= size(v, 1)
 
    ref_data = ref_data';
    
    if size(ref_data, 1) ~= size(v, 1)
        
        display('Reference data must be the same size as simulation data.')
        
        return
        
    end
    
    
end

if center_flag
    
    v = v - mean(v);
    
    ref_data = v - mean(v);
    
end

if norm_flag >= 0
    
    v = v/norm(v, norm_flag);
    
    ref_data = ref_data/norm(ref_data, norm_flag);
    
end

if distance_flag >= 0
   
    diff = norm(v - ref_data, distance_flag);
    
end

results = struct('diff', diff);