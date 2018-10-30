% Example:
%   Define time series:
      y = [102,96.8,97,92.5,95,93,99.4,99.8,105.5];
      t = [97,130,167.5,187.5,218,247.5,285,315,337.5]/360;
%   Define cycle length and alpha:
      w = 2*pi;
      alpha = .05;
%   Run Code:
        t = [t t t t t t];
        y = [y y y y y y];
        figure; plot(t,y)
      cosinor(t,y,w,alpha)
      [cMESOR,cAMPL,cPH,cP] = cosinor2(t,y,1.0);