%% vorticity at outer eyewall boundary
function vort = S_vort_b2(s_i,p) 

    f = p(10);
    beta = p(24);

    vort = f + (1-beta)*(S_v_b2(s_i,p)/S_r_b2(s_i,p));
    end