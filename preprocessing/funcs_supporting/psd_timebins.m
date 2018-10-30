
function [f_bins F_bins t_bins x_bins tcent_bins] = psd_timebins (t,x,timebins)

    if (~exist('fname','var')); fname='blah.bin'; end
    if (~exist('timebins','var')); timebins=2; end % Duration of FFT bins to extract (def 1 hour)
    if (~exist('len','var')); len=Inf; end % Total length of data to read from specified file
    
    plot_on = 0;
%     fs = 12207;
%     if is_downsampled; 
%         fs = fs/12;   % We are using downsampled data
%         len = round(len/12);
%     end


    fs = 1/(mode(diff(t)));
    dt = 1/fs;

    length_bin = round(timebins / dt);
    if (length_bin > length(x))
        fprintf('Timebin size is longer than dataset. Decreasing bin size\n');
        length_bin = length(x);
    end

    
    curr_index = 1;
    nbins = 0;
    
    while (curr_index + length_bin - 1 <= length(x))
        x_bin = x(curr_index:curr_index+length_bin-1);
        x_bin = x_bin - mean(x_bin);
%         t_bin = t(curr_index:curr_index+length_bin-1);
        t_bin = (0:length(x_bin)-1)*dt;

%         [f_bin F_bin] = daveFFT(t_bin,x_bin,1);
        [f_bin F_bin]  = dave_binoverlap_FFT(t_bin, x_bin, 2); % 10 second bin size.
%         [f_bin F_bin] = dave_welch2_FFT(t_bin,x_bin,10);    
        
        N = length(f_bin);
        %T = dt*length(t_bin); %This is wrong because the actual "bin size" is 10 or length of t_bin
        %df = 1/T;
        df = f_bin(2) - f_bin(1);
        T = 1/df;
        dw = 2*pi*df;

        f_bin = f_bin(1:round(N/2)); F_bin = F_bin(1:round(N/2));  % Take only positive half
        F_bin = abs(F_bin).^2 * T / (2*pi)^2 *2;   % Get into mV^2 / Hz (multiply x2 to compensate for -ve half)
        if plot_on;
            
            subplot(211); hold on; plot(t_bin/3600, x_bin,'r'); xlabel('time (h)');
            subplot(212); loglog(f_bin, F_bin,'m'); hold on; 
%             pow_time = sum(abs(x_bin).^2)/length(x_bin)
%             pow_freq = sum(F_bin) * df
%             differen = pow_time - pow_freq

        end;
        nbins = nbins + 1;
        tcent_bins(nbins) = (curr_index + length_bin/2) * dt;
        f_bins{nbins} = f_bin;
        F_bins{nbins} = F_bin;
        t_bins{nbins} = t_bin;
        x_bins{nbins} = x_bin;
        curr_index = curr_index + length_bin;
        
    end
    
    fprintf (['Number of timebins ' num2str(nbins) '.\n']);
    

end