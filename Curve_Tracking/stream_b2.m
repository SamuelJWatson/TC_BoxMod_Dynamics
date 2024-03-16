%% mass stream function at outer eyewall boundary
function stream = stream_b2(s_i,p) 
    
    rho_b = 1.1; % (kg/m^3) Mean boundary layer density
    C_D = 0.003; % Transfer coefficient for momentum

    stream = 2.*pi.*r_b2(s_i,p).*rho_b.*C_D.*(abs(v_b2(s_i,p)).*v_b2(s_i,p)./vort_b2(s_i,p));
end