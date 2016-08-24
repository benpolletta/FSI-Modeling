
function s3 = getAperiodicPulseFast(freq,width,shift,T,dt,onset,offset,ap_pulse_num,ap_pulse_delay,Npop,kernel_type,width2_rise,plot_demo_on)

% Comment this out because it gets confusing.
% if nargin < 5
%     onset = 0;
%     offset = Inf;
% end
% if nargin < 7
%     ap_pulse_num = 0;
%     ap_pulse_delay = 0;
% end
% if nargin < 10
%     kernel_type = 1;  % Options for kernel type are:
%                                 % 1-Gaussian
%                                 % 2-Double exponential - width is decay time; width2 is rise
% end
% if nargin < 11
%     width2_rise = 0.5; % 0.5 ms by default
% end

% plot_demo_on = 1;  % Plot if desired

% Build train of delta functions, spaced with frequency "freq"
t=(0:dt:T)';                            % Generate times vector
s = zeros(size(t));
pulse_period=1000/freq;
s((1+round(shift/dt)):round(pulse_period/dt):end) = 1;    % Add deltas, allowing for shift

% Set aperiodic pulse
if ap_pulse_num > 0
    ap_ind_orig = 1+round(shift/dt)+round(pulse_period/dt)*(ap_pulse_num-1);    % Index of the aperiodic pulse in the time series.
    ap_ind_new = ap_ind_orig+round(ap_pulse_delay/dt);          % Index of where it should appear after the delay.
    if ap_ind_new > length(s) || ap_ind_orig > length(s); error('Aperiodic spike placement would be outside of simulation time.'); end
    s(ap_ind_orig) = 0;                 % Delete the original pulse
    s(ap_ind_new) = 1;                  % Create pulse at the delayed location
end

% Remove anything outside of onset to offset
s(t<onset | t>offset) = 0;


% Build kernel
if kernel_type == 0
        % Build kernel time series
        kernel_length=8*width;                   
        t2a = [0:-dt:-kernel_length];
        t2a = zeros(size(t2a));                  % This applies the heaviside function; no negative times. 
        t2b = [0:dt:kernel_length-dt];
        t2 = [fliplr(t2a(2:end)), 0, t2b(2:end)];
        
        % Build kernel
        kernel = (exp(-t2/width) - exp(-t2/width2_rise)) * (width*width2_rise)/(width-width2_rise);     % This normalization constant is wrong
        kernel = kernel / max(kernel);          % Do normalization manually for now.
elseif kernel_type == -1
        % Build kernel time series
        kernel_length=4*width;                      % Length of kernel time series
        t2a = [0:-dt:-kernel_length];
        t2b = [0:dt:kernel_length-dt];
        t2 = [fliplr(t2a(2:end)), 0, t2b(2:end)];   % This affords us a bit more control over the time values, ensuring it is centered at zero.
        %t2 = -kernel_length:dt:kernel_length-dt;    % Generate time vector
        
        % Build kernel
        kernel = exp(-t2.^2/width^2);      % Build kernel. Peaks at 1.0.
elseif kernel_type > 0
        % Build kernel time series
        kernel_length=4*width;                      % Length of kernel time series
        t2a = [0:-dt:-kernel_length];
        t2b = [0:dt:kernel_length-dt];
        t2 = [fliplr(t2a(2:end)), 0, t2b(2:end)];   % This affords us a bit more control over the time values, ensuring it is centered at zero.
        %t2 = -kernel_length:dt:kernel_length-dt;    % Generate time vector
        
        % Build kernel
        kernel = double(zeros(size(t2)) + abs(t2) < ((width - width/kernel_type)/2 + 1));      % Build smoothing kernel.
        t3a = [0:-dt:-kernel_length];
        t3b = [0:dt:kernel_length-dt];
        t3 = [fliplr(t3a(2:end)), 0, t3b(2:end)];   % This affords us a bit more control over the time values, ensuring it is centered at zero.
        %t2 = -kernel_length:dt:kernel_length-dt;    % Generate time vector
        
        % Build kernel
        smooth_kernel = exp(-t3.^2/(width/kernel_type)^2);      % Build kernel. Peaks at 1.0.
        kernel = conv(kernel, smooth_kernel, 'same');
        
else
        % For debugging; should not reach this!!
        error ('kernel_type should be either 1 (double exponential) or n (Gaussian^n).');
end

kernel = kernel./max(kernel);

kernel = kernel(:);

Nnegatives = length(t2a)-1;
Npositives = length(t2b)-1;
kernel_center = Nnegatives + 1;


% Instead of doing convolution, cycle through the pulse train and find all
% the "1's". Then add the kernel centered at these 1's
N = length(s);
ind = find(s > 0.9);    % Find deltas, record indices
s2=zeros(length(s),1);       % ReSet all deltas to zero; elongate array
for i = 1:length(ind)
    icurr = ind(i);
    istart = max(1,icurr - Nnegatives);     % Dont allow it to be outside array bounds
    istop = min(N,icurr + Npositives);
    
    s2(istart:istop) = s2(istart:istop) + kernel([istart-icurr+kernel_center]:[istop-icurr+kernel_center]);
    
end

% s3=s2;
% 
% % Convolve kernel with deltas
% s2=conv(s,kernel(:));
N2=length(s2); N=length(s);
starting = round((N2-N)/2);
s3=s2(1+starting:N+starting);       % Each edge we're dropping should be half the difference in the length of the two vectors.
%s2=wkeep(s2,length(s),'c');        % wkeep doesn't work with compiled code

if plot_demo_on                 % Plot if desired

    figure; 
    subplot(311); plot(t,s);
    plot_t2=[1:length(s2)]*dt;
    plot_t3=plot_t2(1+starting:N+starting);
    axis tight
    legend('Delta train');
    subplot(312); plot(t2,kernel);
    axis tight
    legend('Kernel');
    subplot(313); plot(plot_t2,s2);
    hold on; plot(plot_t3,s3);
    axis tight
    legend('Original convolution','Keep only center');
    
end


s3 = repmat(s3(:),[1,Npop]);
    
end
