%%Low order TC model 
% clc
% clear

%% Model parameters
g = 9.806; % (m/s^2)
r_a = 420*1000; % (m) Outer radius where p_s = p_ref,s
r_ba = 420*1000; % (m) Outer radius of the ambient region
tau_E = 48; % (h) Timescale for diabatic cooling 
tau_C = 4; % (h) Timescale for convective exchange 
C_H = 0.003; % Transfer coefficient for enthalpy
C_D = 0.003; % Transfer coefficient for momentum
H = 13.5*1000; % (m) Tropopause height - bondary height
H_b = 1.5*1000; % (m) Boundary height
f = 0.00005; % (s^-1) Coriolis parameter
kap = 3; % Eyewall entropy profile parameter
del = 0.25; % Entrainment parameter
R_1 = 90*1000; % (m) Inner potential radius of eyewall
R_2 = 180*1000; % (m) Outer potential radius of eyewall
delR = 30*1000; % (m) Distance from eyewall to outer region
rho = 0.45; % (kg/m^3) Mean density
rho_b = 1.1; % (kg/m^3) Mean boundary layer density
T_t = 203.15; % (K) Tropopause temperature
T_s = 26.8 + 273.15; % (K) Sea surface temperature 
h_a = 0.45; % (%) Relative humidity, ambient region
p_a = 500;
h_refb = 0.80; % (%) Relative humidity, boundary layer
p_ref = 1000;
beta = 0.8; % Tangential wind profile parameter
gamma = (T_s-T_t)/H; % (K/m) temperature lapse rate
M = pi*rho*H*R_2^2; % (kg) total mass contained in angular momentum surface formed by R_2 
M_e = pi*rho*H*R_1^2; % (kg) mass of eye
M_i = M - M_e; % eyewall mass    
L_v = 2264; % (J/g)
R_d = 287; % (J/kgK)
T_a = T_s*(p_a/p_ref)^(R_d*gamma/g); % (K) Temperature ambient region
q_v = SHsat(T_s-273.15); % (g/kg) specific humudity at saturation, need to work out what this should be
q_vas = 1.7.*SHsat(T_a-273.15); % (g/kg) specific humidity  at saturation for ambient region at pressure level p_a, need to work out what this is
q_va = q_vas*h_a; % (g/kg) specific humidity for ambient region at pressure level p_a, need to work out what this is
q_vref = q_v*h_refb; %q_v/2; % (g/kg) specific humudity at p_ref? need to work out what this should be
T_ref = T_s; % (K) reference temperature, what should this be for ambient region?
c_p = 1005; % (J/kg*K)
s_a = L_v*(q_va/T_a - q_vref/T_ref) - R_d*log(p_a/p_ref) + c_p*log(T_a/T_ref); % (J/kg*K) entropy ambient region
s_as = L_v*(q_vas/T_a - q_vref/T_ref) - R_d*log(p_a/p_ref) + c_p*log(T_a/T_ref); % (J/kg*K) saturation entropy ambient region
s_oa0 = L_v*((q_v-q_vref)/T_s); % (J/kg*K) sea surface entropy at surface as R -> inf

% parameter vector
p = [g; %1
    r_a; %2
    r_ba; %3
    tau_E; %4
    tau_C; %5
    C_H; %6
    C_D; %7
    H; %8
    H_b; %9
    f; %10
    kap; %11
    del; %12
    R_1; %13
    R_2; %14
    delR; %15
    rho; %16
    rho_b; %17
    T_t; %18
    T_s; %19 
    h_a; %20
    p_a; %21
    h_refb; %22
    p_ref; %23
    beta; %24
    gamma; %25
    M; %26
    M_e; %27
    M_i; %28    
    L_v; %29
    R_d; %30
    T_a; %31
    q_v; %32
    q_vas; %33
    q_va; %34
    q_vref; %35
    T_ref; %36
    c_p; %37
    s_a; %38
    s_as; %39
    s_oa0]; %40

%% Run Model
%s0 = [20,20,20];
%s0 = [11.592234267,13.2504999935,-2.6963];
%s0 = [2,2,2];
s0 = [-0;-0;-0];
tspan = [0,40];
h = 0.001;

%[t,s] = ode45(@(t,s)TC2(t,s,p),tspan,s0);
%s = s';

% s_i0 = fsolve(@(s_i)v_b2(s_i,p),[0]);
% s_i0
% s0 = [s_i0;0;0];

% s0 = fsolve(@(s)TC2equi(s,p),[s_i0;0;0]);
% s0

[t,s] = RK4(@(t,s) TC2(t,s,p),tspan,h,s0);
s(:,end)

% convert to tangential wind
windy = NaN(1,length(s(1,:)));
for i = 1:length(windy)
   windy(i) = v_b2(s(1,i),p); 
end

%% Plot
seatemp = T_s - 273.15;

figure(1)
clf
til = tiledlayout(2,1);

% plot entropy
ax1 = nexttile;
hold on;
plot(t,s(1,:),'Color',[0 1 0.502],'DisplayName','Eyewall',LineWidth = 5)
plot(t,s(2,:),'Color',[1 0.3 0.3],'DisplayName','Eyewall BL',LineWidth = 5)
plot(t,s(3,:),'Color',[0.2 0.6 1],'DisplayName','Ambient BL',LineWidth = 5) 
ylabel('Entropy (Jkg^{-1}K^{-1})') 
legend()
xticklabels(ax1,{})
ax1 = gca; 
ax1.FontSize = 16;  
% title(sprintf('Entropy and Tangential Wind Evolution for SST = %.f deg C',seatemp))

% plot tangential wind
ax2 = nexttile;
plot(t,windy,'Color',[1 0.2 0.6],LineWidth = 5)
xlabel('Time (hrs)') 
ylabel({'Tangential'; 'wind (ms^{-1})'})  
ax2 = gca; 
ax2.FontSize = 16; 

til.TileSpacing = 'compact';

