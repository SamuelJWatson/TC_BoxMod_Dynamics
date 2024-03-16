%% vorticity at outer eyewall boundary
function vort = vort_b2(s_i,p) 

    f = 0.00005; % (s^-1) Coriolis parameter ****
    beta = p(4,:);

    vort = f+(1-beta).*(v_b2(s_i,p)./r_b2(s_i,p));
    end