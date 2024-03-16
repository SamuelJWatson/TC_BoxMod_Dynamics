%% radius inner eyewall boundary
function r = S_r_b1(s_i,p) 

    gamma = p(25);
    f = p(10);
    R_1 = p(13);
    s_as = p(39);
    delR = p(15);
    R_2 = p(14);
    kap = p(11);
    M_e = p(27);
    rho = p(16);
    H = p(8);

    G_1 = (2*gamma/(f^2*R_1^3))*((s_as-s_i)/delR)*(R_1/R_2)^(kap-1);
        
    r = sqrt((exp(G_1*M_e/(pi*rho))-1)/(G_1*H));
    end