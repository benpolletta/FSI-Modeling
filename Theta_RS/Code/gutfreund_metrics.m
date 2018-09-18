function results = gutfreund_metrics(data, varargin)

if ~isempty(varargin)
    
    if strcmp(varargin{1}, 'variable')
        
        variable = varargin{2};
        
    end
    
end

v = getfield(data, variable);

t = data.time;

rmp = mean(v(t >= 100 & t <= 500));

t_lims = [3001 6000];

v_injected = v(t >= t_lims(1) & t <= t_lims(2));

[v_hat, f] = pmtm(detrend(v_injected), 2, [], 1000*length(t)/t(end), 'eigen');

gauss_kernel = normpdf(-2:.2:2, 0, .5);

gauss_kernel = gauss_kernel/sum(gauss_kernel);

v_hat_smoothed = conv(v_hat, gauss_kernel, 'same');

% v_hat_smoothed = v_hat;

max_freq = f(v_hat_smoothed == max(v_hat_smoothed));

% sto_index = t >= 1001 & t <= 5000;

low = min(v_injected); high = max(v_injected); % low = min(v_injected(sto_index)); high = max(v_injected(sto_index));

rmp_sto = mean(v_injected); % mean(v_injected(sto_index));

amp = high - low;

results = struct('rmp', rmp, 'max_freq', max_freq, 'low', low, 'high', high, 'rmp_sto', rmp_sto, 'amp', amp);