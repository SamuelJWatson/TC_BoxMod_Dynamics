%% mass stream function at outer eyewall boundary
function stream = stream_b2(s_i,p) 
    
    rho_b = p(17);
    C_D = p(7);

    stream = 2*pi*r_b2(s_i,p)*rho_b*C_D*(abs(v_b2(s_i,p))*v_b2(s_i,p)/vort_b2(s_i,p));
end