%% tangential wind speed outer eyewall boundary
function v = v_b2(s_i,p) 

    f = p(10);
    R_2 = p(14);

    v = (f/2)*((R_2^2 - r_b2(s_i,p)^2)/r_b2(s_i,p));
end