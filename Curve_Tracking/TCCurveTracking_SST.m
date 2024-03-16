%% TC model curve tracking (SST)
% Tracks equilibria and bifurcations of TC model through parameter 
% variation using COCO (avaliable at https://sourceforge.net/projects
% /cocotools/).
clc
clear

%% continuation parameters
SST = 28 + 273.15;     % sea surface temp (k) for equilibria
minSST = 16 + 273.15;  % sea surface temp (k) min
maxSST = 32 + 273.15;  % sea surface temp (k) max
beta = 0.8;            % exponent of radial decline
minbeta = 0;           % min beta
maxbeta = 1;           % max beta
h_a = 0.45;            % ambient region relative humidity (%)
h_refb = 0.8;          % boundary layer relative humidity (%)

par = [SST; h_a; h_refb; beta];   

%% initial condition
seq1 = fsolve(@(s)TC2equi_CT(s,par),[0;0;0]); 

%% Bifurcation analysis 
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.02,'h_max',0.1,'h_min',0.001,'ItMX',4000);
bd1 = coco(prob,'eq1','ode','isol','ep',... 
     @(s,par) TC2_coco(s,par),seq1,{'SST' 'h_a' 'h_refb' 'beta'},...
     [SST h_a h_refb beta],...  
     'SST',[minSST maxSST]);
 
 
%% Plot bif diagram
figure(1); clf; hold on;
[SN1,HB1,BP1] = bifdiag1D_SST(bd1,1,'SST');
set(gca,'FontSize',16)
% title('Bifurcation diagram: SST vs Maximum Tangential Wind');
% subtitle(sprintf('h_a = %.2f, h_{ref,b} = %.2f, \\beta = %.3f', h_a, h_refb, beta));
xlabel('SST ^oC');
ylabel('Maximum Tangential Wind (ms^{-1})'); 
xlim([minSST-273.15,maxSST-273.15]) 

%% Continuing fold points in SST and beta
prob = coco_prob();
prob = coco_set(prob,'cont','h0',0.02,'h_max',0.1,'h_min',0.001,'ItMX',2000);

bd_sn31 = coco(prob,'sn31','ode','SN','SN',... 
     'eq1',SN1(1),...  
     {'SST' 'beta'},{[minSST maxSST] [minbeta maxbeta]});
bd_sn32 = coco(prob,'sn32','ode','SN','SN',... 
     'eq1',SN1(2),...  
     {'SST' 'beta'},{[minSST maxSST] [minbeta maxbeta]});

%% Plot 2 param bifurcation continuation for SST and beta
% colours
br = [160 82 45]./255;
dg = [77 149 66]./225;
dp = [75 0  130]./255;
gy = [150 150  150]./255;

SST_sn31 = coco_bd_col(bd_sn31,'SST') - 273.15.*ones(1,length(coco_bd_col(bd_sn31,'SST')));
beta_sn31 = coco_bd_col(bd_sn31,'beta');
SST_sn32 = coco_bd_col(bd_sn32,'SST') - 273.15.*ones(1,length(coco_bd_col(bd_sn32,'SST')));
beta_sn32 = coco_bd_col(bd_sn32,'beta');

figure(2); clf; hold on;
plot(SST_sn31,beta_sn31,'color',dg,'Linewidth',3);
plot(SST_sn32,beta_sn32,'color',dp,'Linewidth',3);
set(gca,'FontSize',16)
title('Saddle-node continuation in \beta and sea surface temperature');
subtitle(sprintf('h_a = %.2f,h_{ref,b} = %.2f', h_a, h_refb));
xlabel('SST ^oC');
ylabel('Exponent of radial decline, \beta'); 