
function s = getPulseFast(t_on,width,T,dt,Npop,kernel_type,plot_demo_on)

% plot_demo_on = 1;  % Plot if desired

% Generate vector of times.
t = (0:dt:T)';         

if ~kernel_type

    s = double((t >= t_on) & (t <= t_on)); s0 = s;
    
elseif (floor(kernel_type) == kernel_type) && isscalar(kernel_type)
    
    s0 = double((t >= t_on + width/(2*kernel_type)) & (t <= t_on + width - width/(2*kernel_type)));
    
    % Build kernel.
    kernel_length = 4*width;                      % Length of kernel time series
    t2a = 0:-dt:-kernel_length;
    t2b = 0:dt:kernel_length-dt;
    t2 = [fliplr(t2a(2:end)), 0, t2b(2:end)];   % This affords us a bit more control over the time values, ensuring it is centered at zero.
    
    smooth_kernel = exp(-t2.^2/(width/kernel_type)^2);      % Build kernel. Peaks at 1.0.
    
    s = conv(s0, smooth_kernel, 'same');
    
else
    
    % For debugging; should not reach this!!
    error ('kernel_type should be an scalar natural number.');
    
end

if plot_demo_on                 % Plot if desired

    figure; 
    subplot(311); plot(t, s0);
    axis tight
    legend('Delta train');
    
    subplot(312); plot(t2, smooth_kernel);
    axis tight
    legend('Kernel');
    
    subplot(313); plot(t, s);
    axis tight
    legend('Pulse');
    
end


s = repmat(s(:),[1,Npop]);
    
end
