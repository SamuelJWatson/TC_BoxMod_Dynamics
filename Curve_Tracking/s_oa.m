%% sea surface entropy under ambient region
function s_ocean_a = s_oa(s_i,p)
    
    L_v = 2264; % (J/g)
    
    T_s = p(1,:);
    h_refb = p(3,:); %0.8; % (%) Relative humidity, boundary layer
    
    q_v = SHsat(T_s-273.15); % (g/kg) specific humudity at saturation, need to work out what this should be
    q_vref = q_v.*h_refb; %q_v/2; % (g/kg) specific humudity at p_ref? need to work out what this should be
 
    s_oa0 = L_v.*((q_v-q_vref)./T_s); % (J/kg*K) sea surface entropy at surface as R -> inf


    s_ocean_a = (s_oi(s_i,p)-s_oa0)./2;
end
