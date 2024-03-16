function dsdt = TC2(t,s,p)
    
    M_i = p(28);
    s_a = p(38);
    s_as = p(39);
    tau_E = p(4);
    C_H = p(6);
    H_b = p(9);
    del = p(12);
    tau_C = p(5);


    dsdt = zeros(3,1); %initialise solution array
    
    %% system of ODEs
    
    % change in eyewall saturation entropy
    dsdt(1) = 40.*(3600*(stream_b2(s(1),p)*(s(2)-s(1))/M_i) + (s_as-s(1))/tau_E);
    
    % change in eyewall boundary layer entropy
    dsdt(2) = 40.*(3600*(stream_b2(s(1),p)*(s(3)-s(2))/M_bi(s(1),p) + (C_H/(2*H_b))*(abs(v_b1(s(1),p))+abs(v_b2(s(1),p)))*(s_oi(s(1),p)-s(2))));
    
    % change in ambient region boundary layer entropy
    dsdt(3) = 40.*(3600*(stream_b2(s(1),p)*(del*s_a-s(3))/M_ba(s(1),p) + (C_H/(2*H_b))*abs(v_b2(s(1),p))*(s_oa(s(1),p)-s(3))) + (s_a-s(3))/tau_C); 
end