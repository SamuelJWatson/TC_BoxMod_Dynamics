%% TC model curve tracking (beta)
% Tracks equilibria and bifurcations of TC model through parameter 
% variation using COCO (avaliable at https://sourceforge.net/projects
% /cocotools/). 
clc
clear

%% continuation parameters
SST = 26.725 + 273.15;   % sea surface temp (k) for equilibria
minSST = 0 + 273.15;     % sea surface temp (k) min
maxSST = 32 + 273.15;    % sea surface temp (k) max
beta = 0.96;             % exponent of radial decline
beta2 = 0.7;
minbeta = 0.5;           % min beta
maxbeta = 0.9999;        % max beta
h_a = 0.45;              % ambient region relative humidity (%)
h_refb = 0.8;            % boundary layer relative humidity (%)


par = [SST; h_a; h_refb; beta];   

%% initial conditions
% upper eq
seq1 = fsolve(@(s)TC2equi_CT(s,par),[0;0;0]); 
seq1

% lower eq
par(4) = beta2;
seq2 = fsolve(@(s)TC2equi_CT(s,par),[0;0;0]); 
seq2

% zero wind state
s_i0 = fsolve(@(s_i)v_b2(s_i,par),[0]);
s_i0

tspan = [0,500];
h = 0.01;

s0 = [s_i0,0,0];

[t,s] = RK4(@(t,s) TC2_CT(t,s,par),tspan,h,s0);

seq3 = fsolve(@(s)TC2equi_CT(s,par),s(:,end));
seq3

%% Bifurcation analysis 
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.005,'h_max',0.01,'h_min',0.001,'ItMX',4000);

% upper eq
bd1 = coco(prob,'eq1','ode','isol','ep',... 
     @(s,par) TC2_coco(s,par),seq1,{'SST' 'h_a' 'h_refb' 'beta'},...
     [SST h_a h_refb beta],...  
     'beta',[minbeta maxbeta]); 

% lower eq
bd2 = coco(prob,'eq2','ode','isol','ep',... 
     @(s,par) TC2_coco(s,par),seq2,{'SST' 'h_a' 'h_refb' 'beta'},...
     [SST h_a h_refb beta2],...  
     'beta',[minbeta maxbeta]);
 
% zero wind eq
bd3 = coco(prob,'eq3','ode','isol','ep',... 
     @(s,par) TC2_coco(s,par),seq3,{'SST' 'h_a' 'h_refb' 'beta'},...
     [SST h_a h_refb beta2],...  
     'beta',[minbeta maxbeta]);
 
%% Plot bif diagram
figure(1); %clf; 
hold on;
[SN1,HB1,BP1] = bifdiag1D(bd1,1,'beta');
[SN2,HB2,BP2] = bifdiag1D(bd2,1,'beta');
[SN3,HB3,BP3] = bifdiag1D(bd3,1,'beta');
set(gca,'FontSize',16);
% title('Bifurcation diagram: beta vs Maximum Tangential Wind');
% subtitle(sprintf('SST = %.2f, h_a = %.2f, h_{ref,b} = %.2f',SST-273.15, h_a, h_refb));
% subtitle(sprintf('SST = %.2f ^oC',SST-273.15));
xlabel('\beta');
ylabel('Maximum Tangential Wind (ms^{-1})');  
xlim([0.5 1]);

%% Continuing fold points in beta and SST
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.02,'h_max',0.1,'h_min',0.001,'ItMX',2000);

bd_sn1 = coco(prob,'sn1','ode','SN','SN',... 
     'eq1',SN1(1),...  
     {'beta' 'SST'},{[minbeta maxbeta] [minSST maxSST]});
bd_sn2 = coco(prob,'sn2','ode','SN','SN',... 
     'eq2',SN2(2),...  
     {'beta' 'SST'},{[minbeta maxbeta] [minSST maxSST]});

%% Plot 2 param bifurcation continuation
% colours
br = [160 82 45]./255;
dg = [77 149 66]./225;
dp = [75 0  130]./255;
gy = [150 150  150]./255;

beta_sn1 = coco_bd_col(bd_sn1,'beta');
SST_sn1 = coco_bd_col(bd_sn1,'SST') - 273.15.*ones(1,length(coco_bd_col(bd_sn1,'SST')));
beta_sn2 = coco_bd_col(bd_sn2,'beta');
SST_sn2 = coco_bd_col(bd_sn2,'SST') - 273.15.*ones(1,length(coco_bd_col(bd_sn2,'SST')));

figure(2); clf; hold on;
plot(beta_sn1,SST_sn1,'color',dg,'Linewidth',3);
plot(beta_sn2,SST_sn2,'color',dp,'Linewidth',3);
set(gca,'FontSize',16)
title('Saddle-node continuation in \beta and SST');
subtitle(sprintf('h_a = %.2f, h_{ref,b} = %.2f',h_a, h_refb));
xlabel('\beta');
ylabel('SST ^oC'); 