
% Function used to construct FileInfo.mat files based on info in original seizure logs.
function build_fileInfo (ratN)

    if ratN == 4
        
        load R004FileInfo_orig.mat
        
        Sz = [];
        Sz = [Sz; 99 calc_offset(1,35,0)];  % Note = Previous versions of FileInfo do not take into account seconds, and therefore we do not include them here.
        Sz = [Sz; 192 calc_offset(0,59,0)];
        Sz = [Sz; 209 calc_offset(3,5,0)];
        Sz = [Sz; 216 calc_offset(0,16,0)];
        Sz = [Sz; 236 calc_offset(1,37,0)];
        
        save ('R004FileInfo.mat','Sz','UB','LB')
        
    elseif ratN == 1
        
        Sz = [];
        LB = [];
        UB = [];
        Sz = [Sz; 109 calc_offset(7,19,0)];
        Sz = [Sz; 114 calc_offset(7,14,0)];
        Sz = [Sz; 115 calc_offset(7,15,0)];
        Sz = [Sz; 126 calc_offset(4,37,0)];
        Sz = [Sz; 127 calc_offset(4,18,0)];
        Sz = [Sz; 168 calc_offset(1,8,0)];
        
        save ('R001FileInfo.mat','Sz','UB','LB')
    end

end

function t_sec = calc_offset(h,min,sec)
    t_sec = h*3600+min*60+sec*1;
end