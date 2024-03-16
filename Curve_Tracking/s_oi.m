%% sea surface entropy under eyewall, s_oi = so2 (sea surface entropy at outer eyewall boundary)
function s_ocean_i = s_oi(s_i,p) 

    L_v = 2264; % (J/g)
    f = 0.00005; % (s^-1) Coriolis parameter ****
    r_a = 420*1000; % (m) Outer radius where p_s = p_ref,s
    
    T_s = p(1,:);
    h_refb = p(3,:);
    beta = p(4,:);
    
    q_v = SHsat(T_s-273.15); % (g/kg) specific humudity at saturation, need to work out what this should be
    q_vref = q_v.*h_refb; %q_v/2; % (g/kg) specific humudity at p_ref? need to work out what this should be
 

    s_ocean_i = L_v.*((q_v-q_vref)./T_s) + (v_b2(s_i,p).^2./(T_s.*2.*beta)).*(1-(r_b2(s_i,p)./r_a).^(2.*beta)) - ...
        (f.*v_b2(s_i,p).*r_b2(s_i,p)./(T_s.*(1-beta))).*(1-(r_a./r_b2(s_i,p)).^(1-beta));
    end