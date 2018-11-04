function TM = buildTimeMatrix(Nblocks,tstart,tjump,tskip)
% Builds matrix TM
% Each row is a time range for plotting

if nargin < 4
    tskip = 0;
end

TM = [tstart+tjump*0+tskip*0:(tjump+tskip):tstart+(tjump+tskip)*(Nblocks-1); ...
      tstart+tjump*1:(tjump+tskip):tstart+(tjump+tskip)*Nblocks]';

end