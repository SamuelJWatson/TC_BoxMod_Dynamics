%% mass of eyewall boundary layer
function mass_bi = M_bi(s_i,p) 

    rho_b = 1.1; % (kg/m^3) Mean boundary layer density
    H_b = 1.5*1000; % (m) Boundary height   

    mass_bi = pi.*rho_b.*(r_b2(s_i,p).^2-r_b1(s_i,p).^2).*H_b;
    end