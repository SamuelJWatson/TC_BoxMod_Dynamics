%% mass of ambient boundary layer
function mass_ba = S_M_ba(s_i,p) 

    rho_b = p(17); 
    r_ba = p(3);
    H_b = p(9);

    mass_ba = pi*rho_b*(r_ba^2-S_r_b2(s_i,p)^2)*H_b;
    end