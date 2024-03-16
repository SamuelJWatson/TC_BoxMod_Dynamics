%% tangential wind speed inner eyewall boundary
function v = v_b1(s_i,p) 

    f = 0.00005; % (s^-1) Coriolis parameter ****
    R_1 = 90*1000; % (m) Inner potential radius of eyewall

    v = (f/2).*((R_1.^2-r_b1(s_i,p).^2)./r_b1(s_i,p));
    end