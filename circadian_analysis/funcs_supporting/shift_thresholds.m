
function phi = shift_thresholds(phi,thresh)

    phi = abs(phi/2/pi);    % Take abs value and normalize radians to 1.0 ("days").
    phi(phi < thresh) = phi(phi < thresh) + 1;
    
end