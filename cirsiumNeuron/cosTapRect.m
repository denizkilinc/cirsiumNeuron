function [y, yp] = cosTapRect(t, T, r, K, A)
    % COSTAPRECT cosine tapered pulse signal.

    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------

    % [y yp] = COSTAPRECT(t, T, r, K, A)
    % t... time
    % T... signal period
    % r... ratio of cosine ramp width to pulse width
    % K... ratio of pulse width to signal duration T
    % A... pulse amplitude
    t = (t/T - floor(t/T))/K;
    if t < r
        y = 1/2*(1 + cos((pi/r)*(t-r)));
        yp = -1/2*sin((pi/r)*(t-r))*pi/r;
    elseif t < (1-r)
        y = 1;
        yp = 0;
    elseif t < 1
        y = 1/2*(1 + cos(pi/r*(t-1+r)));
        yp = -1/2*sin(pi/r*(t-1+r))*pi/r;
    else
        y = 0;
        yp = 0;
    end
    y = A*y;
    yp = A*yp/T/K;
end