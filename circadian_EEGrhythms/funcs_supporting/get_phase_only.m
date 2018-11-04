function phi = get_phase_only(t,x,w,thresh,amplitude)

    if nargin < 5
        amplitude = 0;
    end

%     w = 2*pi/1.0;
    alpha = .05;

    sout = cosinor_struct(t,x,w,alpha);
    
    phi = sout.phi;
    
    phi = shift_thresholds(phi,thresh);
%     sout.phi = phi;
%     sout.CI_phi_min = CI_phi_min;
%     sout.CI_phi_max = CI_phi_max;
%     sout.RNE = RNE;
%     sout.p_3a = p_3a;
    if amplitude
        phi = sout.Amp;  %output amplitude instead of phase        
    end
end
