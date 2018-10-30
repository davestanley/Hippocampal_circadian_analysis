

function [Fout] = extract_FFT_band (f,F,freq_band)
    
%     f = struct2matrix(file,'f',[],0,0);
%     F = struct2matrix(file,'F',[],0,0);
%     tabs = struct2matrix(file,'tabs',[],0,0);
    

    freq_band = freq_band(:);
    
    index = find((f(:,1) >= freq_band(1)) & (f(:,1) < freq_band(2)));
    f_width = abs(diff(freq_band));
    
    Fout = sum(F( index,: ),1);
    Fout = Fout / f_width;

end
