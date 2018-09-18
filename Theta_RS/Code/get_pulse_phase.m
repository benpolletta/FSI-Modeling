function results = get_pulse_phase(data, varargin)

if ~isempty(varargin)
    
    if strcmp(varargin{1}, 'variable')
        
        variable = varargin{2};
        
    elseif strcmp(varargin{1}, 'figure_flag')
        
        figure_flag = varargin{2};
        
    end
    
else
    
    variable = 'pop1_v'; figure_flag = 0;
    
end

v = getfield(data, variable);

t = data.time;

%% Getting pulse times, and voltage time series before and after pulse is turned on.

pulse_on = data.pop1_iPulse_PPon;

pulse_width = data.pop1_iPulse_PPwidth;

pulse_off = pulse_on + pulse_width;

v_before = v(t < pulse_on);

v_after = v(t > pulse_off);

%% Getting peak frequency.

[v_hat, f] = pmtm(detrend(v), 2, [], sampling_freq, 'eigen');

gauss_kernel = normpdf(-2:.2:2, 0, .5);

gauss_kernel = gauss_kernel/sum(gauss_kernel);

v_hat_smoothed = conv(v_hat, gauss_kernel, 'same');

peak_freq = f(v_hat_smoothed == max(v_hat_smoothed));

%% Getting bandpassed time series.

sampling_freq = 1000*length(t)/t(end);

v_bandpassed_before = wavelet_spectrogram(v_before, sampling_freq, peak_freq, 7, 0, '');

v_bandpassed_after = wavelet_spectrogram(v_after, sampling_freq, peak_freq, 7, 0, '');

v_phases = [angle(v_bandpassed_before((end - 100))) angle(v_bandpassed_after(100))];

if figure_flag
   
    subplot(2, 1, 1)
    
    ax = plotyy(t, v, t(t < pulse_on), real(v_bandpassed_before));
    
    axis(ax, 'tight')
    
    subplot(2, 1, 2)
    
    ax = plotyy(t, v, t(t < pulse_on), real(v_bandpassed_before));
    
    axis(ax, 'tight')
    
end
    
end