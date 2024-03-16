%% radius inner eyewall boundary
function r = r_b1(s_i,p) 

    f = 0.00005; % (s^-1) Coriolis parameter ****
    R_1 = 90*1000; % (m) Inner potential radius of eyewall
    delR = 30*1000; % (m) Distance from eyewall to outer region
    R_2 = 180*1000; % (m) Outer potential radius of eyewall
    kap = 3; % Eyewall entropy profile parameter
    rho = 0.45; % (kg/m^3) Mean density
    H = 13.5*1000; % (m) Tropopause height - bondary height
    M_e = pi.*rho.*H.*R_1.^2; % (kg) mass of eye
    T_t = 203.15; % (K) Tropopause temperature
    
    L_v = 2264; % (J/g)
    R_d = 287; % (J/kgK)
    c_p = 1005; % (J/kg*K)
    p_a = 500;
    p_ref = 1000;
    g = 9.806; % (m/s^2)
    
    T_s = p(1,:);
    h_refb = p(3,:); %0.8; % (%) Relative humidity, boundary layer
    gamma = (T_s-T_t)./H; % (K/m) temperature lapse rate
    
    T_a = T_s.*(p_a./p_ref).^(R_d.*gamma./g); % (K) Temperature ambient region
    q_v = SHsat(T_s-273.15); % (g/kg) specific humudity at saturation, need to work out what this should be
    q_vas = 1.7.*SHsat(T_a-273.15); % (g/kg) specific humidity  at saturation for ambient region at pressure level p_a, need to work out what this is
    q_vref = q_v.*h_refb; %q_v/2; % (g/kg) specific humudity at p_ref? need to work out what this should be
    T_ref = T_s; % (K) reference temperature, what should this be for ambient region?
    
    s_as = L_v.*(q_vas./T_a - q_vref./T_ref) - R_d.*log(p_a./p_ref) + c_p.*log(T_a./T_ref); % (J/kg*K) saturation entropy ambient region
    
    G_1 = (2.*gamma./(f.^2.*R_1.^3)).*((s_as-s_i)./delR).*(R_1./R_2).^(kap-1);
        
    r = sqrt((exp(G_1.*M_e./(pi.*rho))-1)./(G_1.*H));
    end