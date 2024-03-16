%% mass of eyewall boundary layer
function mass_bi = S_M_bi(s_i,p) 

    rho_b = p(17);
    H_b = p(9);    

    mass_bi = pi*rho_b*(S_r_b2(s_i,p)^2-S_r_b1(s_i,p)^2)*H_b;
    end