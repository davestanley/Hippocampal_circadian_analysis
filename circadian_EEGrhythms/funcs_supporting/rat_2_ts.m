

function [T X] = rat_2_ts (rcell,prepostchronic, smoothmode_ts)

    for j = 1:length(rcell)
        r=rcell{j}; [T{j} X{j}] = extract_ts(r, prepostchronic, smoothmode_ts); 
    end
end