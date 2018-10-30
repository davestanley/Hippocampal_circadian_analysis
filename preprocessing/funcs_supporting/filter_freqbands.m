

function [xarr] = filter_freqbands (t,x,freq_band)

    xarr = zeros(length(t),size(freq_band,1));
    for i = 1:size(freq_band,1)
        xarr(:,i) = qif(t,x,[0 freq_band(i,1); freq_band(i,2) Inf]);
    end
    
end


