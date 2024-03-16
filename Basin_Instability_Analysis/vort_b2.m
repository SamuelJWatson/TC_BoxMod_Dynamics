%% vorticity at outer eyewall boundary
function vort = vort_b2(s_i,p) 

    f = p(10);
    beta = p(24);

    vort = f+(1-beta)*(v_b2(s_i,p)/r_b2(s_i,p));
    end