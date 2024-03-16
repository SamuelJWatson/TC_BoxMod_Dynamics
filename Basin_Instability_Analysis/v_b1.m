%% tangential wind speed inner eyewall boundary
function v = v_b1(s_i,p) 

    f = p(10);
    R_1 = p(13);

    v = (f/2)*((R_1^2-r_b1(s_i,p)^2)/r_b1(s_i,p));
    end