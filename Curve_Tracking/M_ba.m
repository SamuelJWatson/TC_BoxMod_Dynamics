%% mass of ambient boundary layer
function mass_ba = M_ba(s_i,p) 

    r_ba = 420*1000; % (m) Outer radius of the ambient region
    rho_b = 1.1; % (kg/m^3) Mean boundary layer density
    H_b = 1.5*1000; % (m) Boundary height

    mass_ba = pi.*rho_b.*(r_ba.^2-r_b2(s_i,p).^2).*H_b;
    end