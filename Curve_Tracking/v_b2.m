%% tangential wind speed outer eyewall boundary
function v = v_b2(s_i,p) 

    f = 0.00005; % (s^-1) Coriolis parameter ****
    R_2 = 180*1000; % (m) Outer potential radius of eyewall

    v = (f/2).*((R_2.^2 - r_b2(s_i,p).^2)./r_b2(s_i,p));
    %disp(v*1000)
end