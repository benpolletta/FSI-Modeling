function results = phase_metrics(data, varargin, figure_flag)

if ~isempty(varargin)
    
    if strcmp(varargin{1}, 'variable')
        
        variable = varargin{2};
        
    end
    
else
    
    variable = 'pop1_v';
    
end

v = getfield(data, variable);

t = data.time;

input = data.pop1_iPeriodicPulses_Iext;

sampling_freq = 1000*length(t)/t(end);

%% Getting peak frequency.

[v_hat, f] = pmtm(detrend(v), 2, [], sampling_freq, 'eigen');

gauss_kernel = normpdf(-2:.2:2, 0, .5);

gauss_kernel = gauss_kernel/sum(gauss_kernel);

v_hat_smoothed = conv(v_hat, gauss_kernel, 'same');

max_freq = f(v_hat_smoothed == max(v_hat_smoothed));

freqs = [1.5 4.5 max_freq]'; no_cycles = [7 7 7]'; no_freqs = length(freqs);

%% Getting wavelet components.

v_bandpassed(:, 1) = wavelet_spectrogram(-input, sampling_freq, freqs(1), no_cycles(1), 0, '');

v_bandpassed(:, 2:3) = wavelet_spectrogram(v, sampling_freq, freqs(2:3), no_cycles(2:3), 0, '');

v_phase = angle(v_bandpassed);

if figure_flag

    figure, subplot(4, 1, 1), plot(f(f <= 15), v_hat(f <= 15))
    
    t_end = t >= 0; % 3000;
    
    for frequency = 1:3, 
        
        subplot(4, 2, 2 + 2*(frequency - 1) + 1)
        
        plotyy(t(t_end), real(v_bandpassed(t_end, frequency)), t(t_end), v(t_end))
    
        subplot(4, 2, 2 + 2*(frequency - 1) + 2)
        
        plotyy(t(t_end), angle(v_bandpassed(t_end, frequency)), t(t_end), v(t_end))
        
    end
    
end

%% Getting spike times and computing spike phases.

v_spikes = [diff(v > 0) == 1; zeros(size(v, 2))];

v_spike_phases = v_phase(repmat(logical(v_spikes), 1, no_freqs));

if figure_flag
   
    figure
    
    subplot(2, 1, 1)
    
    plotyy(t, v, t, v_spikes)
    
    for i = 1:3
        
        subplot(2, 3, 3 + i)
        
        rose(gca, v_spike_phases(:, i))
    
    end
    
end

%% Computing phase-phase relationships.

no_bins = 18;

bin_centers = 1:no_bins;

bin_centers = (bin_centers - 1)*2*pi/no_bins - pi*(no_bins - 1)/no_bins;

bin_left_endpoints = bin_centers - pi/no_bins;

bin_right_endpoints = bin_centers + pi/no_bins;

f_pairs = nchoosek(1:no_freqs, 2); no_f_pairs = size(f_pairs, 1);

[v_phase_coh, v_phase_angle] = deal(nan(no_bins, no_freqs, no_freqs));

v_phase_phase = nan(no_bins, no_bins, no_f_pairs);

for fp = 1:no_f_pairs
    
    for b1 = 1:no_bins
        
        b1_indicator = v_phase(:, f_pairs(fp, 1)) >= bin_left_endpoints(b1) & v_phase(:, f_pairs(fp, 1)) < bin_right_endpoints(b1);
        
        v_phase_coh(b1, f_pairs(fp, 2), f_pairs(fp, 1)) = circ_r(v_phase(b1_indicator, f_pairs(fp, 2)));
        
        v_phase_angle(b1, f_pairs(fp, 2), f_pairs(fp, 1)) = circ_mean(v_phase(b1_indicator, f_pairs(fp, 2)));
        
        for b2 = 1:no_bins
   
            b2_indicator = v_phase(:, f_pairs(fp, 2)) >= bin_left_endpoints(b2) & v_phase(:, f_pairs(fp, 2)) < bin_right_endpoints(b2);
        
            v_phase_coh(b2, f_pairs(fp, 1), f_pairs(fp, 2)) = circ_r(v_phase(b2_indicator, f_pairs(fp, 1)));
        
            v_phase_angle(b2, f_pairs(fp, 1), f_pairs(fp, 2)) = circ_mean(v_phase(b2_indicator, f_pairs(fp, 1)));
            
            v_phase_phase(b2, b1, fp) = sum(b1_indicator & b2_indicator);
    
        end
        
    end
    
end

results = struct('v_spike_phases', v_spike_phases, 'v_phase_coh', v_phase_coh, 'v_phase_angle', v_phase_angle, 'v_phase_phase', v_phase_phase);