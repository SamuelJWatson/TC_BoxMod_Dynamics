%% basin instability 
% Calculate basin of instability for TC model. 
% Runs model for range of beta an SST with set intial condintion and test
% for convergence to low or high wind equilibria. Then uses a bisection
% method to define boundary line. 
clc
clear

%% TC model
% variation parameters
T_s = 27.1 + 273.15;     % IC SST
SSTmin = 26.4 + 273.15;  % min SST tested
SSTmax = 27.2 + 273.15;  % max SST tested
beta = 0.95;             % IC beta
betamin = 0.55;          % min beta tested
betamax = 1;             % max beta tested
its = 10;               % discretisation

icpoint = [beta,T_s-273.15,100]; % with dummy value for plotting

% Model parameters
g = 9.806; % (m/s^2)
r_a = 420*1000; % (m) Outer radius where p_s = p_ref,s
r_ba = 420*1000; % (m) Outer radius of the ambient region
tau_E = 48; % (h) Timescale for diabatic cooling 
tau_C = 4; % (h) Timescale for convective exchange 
C_H = 0.003; % Transfer coefficient for enthalpy
C_D = 0.003; % Transfer coefficient for momentum
H = 13.5*1000; % (m) Tropopause height - bondary height
H_b = 1.5*1000; % (m) Boundary height
f = 0.00005; % (s^-1) Coriolis parameter ****
kap = 3; % Eyewall entropy profile parameter
del = 0.25; % Entrainment parameter
R_1 = 90*1000; % (m) Inner potential radius of eyewall
R_2 = 180*1000; % (m) Outer potential radius of eyewall
delR = 30*1000; % (m) Distance from eyewall to outer region
rho = 0.45; % (kg/m^3) Mean density
rho_b = 1.1; % (kg/m^3) Mean boundary layer density
T_t = 203.15; % (K) Tropopause temperature
h_a = 0.45; % (%) Relative humidity, ambient region
p_a = 500;
h_refb = 0.80; % (%) Relative humidity, boundary layer
p_ref = 1000;
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

%% Initial condition 
% entropy values of max wind state for SST_0 and beta_0 
% within region of bi-stability.

% find stable solution
s0 = fsolve(@(s)TC2equi(s,p),[0;0;0]);
s0

%% run model to check it is max wind eq
tspan = [0,40];
h = 0.001;

[t,s] = RK4(@(t,s) TC2(t,s,p),tspan,h,s0);

% convert to tangential wind
windy = NaN(1,length(s(1,:)));
for i = 1:length(windy)
   windy(i) = v_b2(s(1,i),p); 
end

% plot
figure(1); clf;
plot(t,windy,'m')
xlabel('Time (hrs)') 
ylabel('Tangential wind (m/s)') 

%% Run with initial condition for range of SST, beta

xxbeta = linspace(betamin, betamax, its);
xxSST = linspace(SSTmin, SSTmax, its);
basinstab = NaN(length(xxSST),length(xxbeta));
windyfinal = NaN(length(xxSST),length(xxbeta));

for i = 1:length(xxbeta)
    p(24) = xxbeta(i); % update beta
    
    sprintf('%d out of %d',i,length(xxbeta))
   for j = 1:length(xxSST)
       
       p(19) = xxSST(j); % update SST (T_s)
       p(25) = (p(19)-p(18))/p(8); % update lapse rate (gamma)
       p(31) = p(19)*(p(21)/p(23))^(p(30)*p(25)/p(1)); % update temperature ambient region (T_a)
       p(32) = SHsat(p(19)-273.15); % update specific humidity at saturation for SST (q_v)
       p(33) = 1.7.*SHsat(p(31)-273.15); % update specific humidity at saturation for ambient region (q_vas)
       p(34) = p(33)*p(20); % update specific humidity for ambient region (q_va)
       p(35) = p(32)*p(22); % update specific humidity for reference (q_vref)
       p(36) = p(19); % update T_ref
       p(38) = p(29)*(p(34)/p(31) - p(35)/p(36)) - p(30)*log(p(21)/p(23)) + p(37)*log(p(31)/p(36)); % update ambient region entropy (s_a)
       p(39) = p(29)*(p(33)/p(31) - p(35)/p(36)) - p(30)*log(p(21)/p(23)) + p(37)*log(p(31)/p(36)); % update ambient region saturation entropy (s_as)
       p(40) = p(29)*((p(32)-p(35))/p(19)); % update sea surface entropy (s_oa0)

       % run model with set IC
       [t,s] = RK4(@(t,s) TC2(t,s,p),tspan,h,s0);
       
       % convert final entropy value to tangential wind
       windyfinal(j,i) = v_b2(s(1,end),p);
       
       % assess whether high of low eq 
       if isnan(windyfinal(j,i))
           basinstab(j,i) = NaN; % model breaks
       elseif windyfinal(j,i) > 15
           basinstab(j,i) = 1; % high wind eq
       else
           basinstab(j,i) = 2; % low wind eq
       end   
   end
end

%% load SST-beta SN continuation data

load('SSTbetaSNCont.mat','SST_sn31','beta_sn31');
vals = 100.*ones(1,length(SST_sn31)); %dummy values

%% plot results

xxctemp = xxSST - 273.15.*ones(1,length(xxSST));

figure(2); clf;
surf(xxbeta,xxctemp,windyfinal);
colorbar;
hold on; 
scatter3(icpoint(1),icpoint(2),icpoint(3),150,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]); %IC point
plot3(beta_sn31,SST_sn31,vals,linewidth = 3); %overlay SN continuation
view(0,90);
ylabel('SST ^oC');
xlabel('\beta');
ylim([26.4,27.2]);
xlim([0.55,1]);
hold off;

%% Bisect method to divide basins of convergence

% make NaNs = 0
for i = 1:size(xxbeta,2)
    for j = 1:size(xxctemp,2)
        if isnan(basinstab(j,i))
            basinstab(j,i) = 0;
        end
    end
end

% Separate basins from collist data matrix
collist = unique(basinstab);
Basin1 = []; % NaNs
Basin2 = []; % low wind state
Basin3 = []; % high wind state

for i = 1:size(xxbeta,2)
    for j = 1:size(xxctemp,2)
        if basinstab(j,i) == collist(1) 
            Basin1 = [Basin1;[xxbeta(i),xxctemp(j)]];
        elseif basinstab(j,i) == collist(2)
            Basin2 = [Basin2;[xxbeta(i),xxctemp(j)]];
        elseif basinstab(j,i) == collist(3)
            Basin3 = [Basin3;[xxbeta(i),xxctemp(j)]];
        else
            print('error - not enough Basins') 
        end
    end
end

% get sizes
Bsize1 = size(Basin1,1);
Bsize2 = size(Basin2,1);
Bsize3 = size(Basin3,1);

% take smallest basin size
bisectsize = min(Bsize2,Bsize3);

A = NaN(bisectsize,2);
B = NaN(bisectsize,2);

% use bisection method to find basin boundary
for i=1:bisectsize
    [A(i,:),B(i,:)] = bisect(Basin2(i,:),Basin3(i,:),Basin2,Basin3,1e-5);
end

%% Sort Boundary data and smooth 
[Axsort, Axind] = sort(A(:,1));
Aysort = A(Axind,2);
Asmooth = smoothdata([Axsort,Aysort],'gaussian',bisectsize/20);

[Bxsort, Bxind] = sort(B(:,1));
Bysort = B(Bxind,2);
Bsmooth = smoothdata([Bxsort,Bysort],'gaussian',bisectsize/20);

%% plot
%colours
blu = [0 128 255]./255;
gre = [0 255 128]./255;
red = [255 51 51]./255;
br = [160 82 45]./255;
dg = [77 149 66]./225;
gra = [160 160 160]./225;

figure(2); clf; hold on;
plot(beta_sn31(1:394),SST_sn31(1:394),'k',linewidth = 4)
plot(Bsmooth(:,1),Bsmooth(:,2),'.','color',blu,'MarkerSize',18)
plot(icpoint(1),icpoint(2),'+','color',blu,'MarkerSize',18,linewidth = 4)
xlabel('\beta') 
ylabel({'SST (^oC)'}) 
ax = gca; 
ax.FontName = 'Calibri';
ax.FontSize = 22; 
box on;
axis([Bsmooth(1,1) Bsmooth(end,1) 26.4 27.2])