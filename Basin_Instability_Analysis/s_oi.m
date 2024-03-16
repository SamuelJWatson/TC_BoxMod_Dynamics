%% sea surface entropy under eyewall, s_oi = so2 (sea surface entropy at outer eyewall boundary)
function s_ocean_i = s_oi(s_i,p) 

    L_v = p(29);
    q_v = p(32);
    q_vref = p(35);
    T_s = p(19);
    beta = p(24);
    r_a = p(2);
    f = p(10);
    

    s_ocean_i = L_v*((q_v-q_vref)/T_s) + (v_b2(s_i,p)^2/(T_s*2*beta))*(1-(r_b2(s_i,p)/r_a)^(2*beta)) - (f*v_b2(s_i,p)*r_b2(s_i,p)/(T_s*(1-beta)))*(1-(r_a/r_b2(s_i,p))^(1-beta));
    end