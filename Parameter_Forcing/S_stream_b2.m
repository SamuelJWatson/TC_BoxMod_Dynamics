%% mass stream function at outer eyewall boundary
function stream = S_stream_b2(s_i,p) 
    
    rho_b = p(17);
    C_D = p(7);

    stream = 2*pi*S_r_b2(s_i,p)*rho_b*C_D*(abs(S_v_b2(s_i,p))*S_v_b2(s_i,p)/S_vort_b2(s_i,p));
end