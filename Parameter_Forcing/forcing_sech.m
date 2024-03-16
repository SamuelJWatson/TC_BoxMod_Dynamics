%% Sech forcing function
function FSech = forcing_sech(p, t)

up = p(1); % decreasing (0) or increasing (1) profile
ramp = p(2); % whether return (0) or ramp (1) profile
rate = p(3); % parameter rate of change
T = p(4); % time at which peak is reached


    if up == 0 % decreasing profile
        if ramp == 0 % return profile
            FSech = -sech(rate*(t-T)) + 1;
        elseif ramp == 1 % ramping profile
            if t <= T 
                FSech = -sech(rate*(t-T)) + 1;
            elseif t > T
                FSech = 0;  % after peak stays constant at max value.
            end        
        end    
    elseif up == 1 % increasing profile
        if ramp == 0 % return profile
            FSech = sech(rate*(t-T));
        elseif ramp == 1 % ramping profile
            if t <= T 
                FSech = sech(rate*(t-T));
            elseif t > T
                FSech = 1;  % after peak stays constant at max value.
            end        
        end
    end
 
end