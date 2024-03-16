%% radius outer eyewall boundary
function r = r_b2(s_i,p) 

    gamma = p(25);
    f = p(10);
    s_as = p(39);
    delR = p(15);
    R_2 = p(14);
    M = p(26);
    rho = p(16);
    H = p(8);

    G_2 = (2*gamma/(f^2*R_2^3))*((s_as-s_i)/delR);
        
    r = sqrt((exp(G_2*M/(pi*rho))-1)/(G_2*H));
    end