

function plot_thetahighlight(t,x,theta_delta_ratio)
    freq_band = [1 5; 5 10; 12 25; 25 50; 50 90; 90 140; 140 200; 200 500]';
    timebin = 2;    % FFT window duration in seconds
    [f F t2 x2 tcent] = psd_timebins (t,x,timebin);
    fout = zeros(size(freq_band,2),length(f));
    Fout = fout;

    for i = 1:length(F)
        [fout(:,i) Fout(:,i)] = psd_freqbins (f{i},F{i}, freq_band); % Bin spectrogram
    end

    [delta] = extract_FFT_band (fout,Fout,freq_band(:,1));
    [theta] = extract_FFT_band (fout,Fout,freq_band(:,2));

    index = theta./delta > theta_delta_ratio;
    hold on; plot(tcent(index),+5.5,'rv','MarkerSize',15,'LineWidth',2)

    tmin = tcent(index)-1.0;
    tmax = tcent(index)+1.0;
    index = logical(zeros(1,length(t)));
    for i = 1:length(tmin)
        index2 = (t >= tmin(i) & t < tmax(i));
        index = ( index |  index2);
        plot(t(index2),x(index2)+4.5,'r')
    end
end