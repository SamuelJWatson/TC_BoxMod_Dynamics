%% Low order TC model with beta and sech forcing 
clc
clear

%% beta sech forcing parameters
% single or double sech forcing profile. To use single profile set time of
% change to greater than model run time.
bsfp = [0.74;    %1 min beta
        0.95;    %2 max beta
        120;     %3 time of change between profiles
        % first forcing
        0;       %4 decreasing(0)/increasing(1) profile
        0;       %5 return(0)/ramp(1) profile
        0.2;     %6 rate
        24;      %7 time of peak forcing
        % second forcing
        0;       %8 decreasing(0)/increasing(1) profile
        0;       %9 return(0)/ramp(1) profile
        1;       %10 rate
        130];    %11 time of peak forcing

%% SST sech forcing parameters
% single sech forcing profile
Ssfp = [26.5 + 273.15; % min SST
        26.9 + 273.15; % max SST
        0;             % decreasing(0)/increasing profile(1)
        1;             % return(0)/ramp profile(1)
        0.3;           % rate
        36];          % time of peak forcing

%% Irma SST Profile
% For use with Irma forcing data
% load('IrmaSSTProfSmooth','Irmat','IrmaSST')

%% Run Model
s0 = [-66.6739;-66.2013;-75.0304]; % initial condition
tspan = [0,60]; 
h = 0.001;           % time step

[t,s,windy,radius] = RK4(@(t,s) S_TC(t,s,bsfp,Irmat,IrmaSST),tspan,h,s0);

%% Get forcing profiles to plot
betaforced = zeros(1,length(t));

for i = 1:length(t)
    betaforced(i) = beta_sech(bsfp, t(i));
end

SSTforced = zeros(1,length(t));

for i = 1:length(t)
    SSTforced(i) = SST_sech(Ssfp, t(i)) - 273.15;
end

% For use with Irma forcing data
% for i = 1:length(t)
%     SSTforced(i) = forcing_IrmaSST(t(i),Irmat,IrmaSST) - 273.15;
% end

rad = radius./(10^4);

%% Plot

% plot entropy
figure(1);
clf;
hold on;
plot(t,s(1,:),'g','DisplayName','Eyewall')
plot(t,s(2,:),'r','DisplayName','Eyewall BL')
plot(t,s(3,:),'b','DisplayName','Ambient BL')
xlabel('Time (hrs)') 
ylabel('Entropy (J/kgK)') 
% title(sprintf('Entropy with \\beta forced by %s %s profile and SST by a %s profile',direction,type,Sdirection), ...
%     sprintf('SST = [%.3f,%.3f] %cC, \\beta = [%.3f,%.3f], rate = %.2f, peak = %d hrs.',Ssfp(1),Ssfp(2),char(176),bsfp(1),bsfp(2),bsfp(5),bsfp(6)))
legend()

% plot tangential wind
figure(2);
clf;
til = tiledlayout(3,1);

ax1 = nexttile;
plot(t,windy,'Color',[1 0.2 0.6],LineWidth = 5)
ylabel({'Tangential'; 'wind (ms^{-1})'}) 
% title('Tangential wind')
ylim([min(windy)-5, max(windy)+5])
xlim([tspan(1),tspan(2)])
xticklabels(ax1,{})
ax1 = gca; 
ax1.FontSize = 16; 

ax2 = nexttile;
plot(t,rad,'Color',[1 0.3 0.3],LineWidth = 5)
ylabel({'Radius'; '\times 10^{4}(m)'}) 
ylim([min(rad)-0.2, max(rad)+0.2])
xlim([tspan(1),tspan(2)])
% title('Outer eyewall raduis')
xticklabels(ax2,{})
ax2 = gca; 
ax2.FontSize = 16; 

ax3 = nexttile;
plot(t,betaforced,'Color',[0.2 0.6 1],LineWidth = 5)
% xlabel('Time (hrs)')
ylabel('\beta') 
ylim([min(betaforced)-0.05, max(betaforced)+0.05])
xlim([tspan(1),tspan(2)])
% title('\beta forcing profile',...
%     sprintf('\\beta = [%.2f,%.3f], rate = %.2f, peak = %d hrs.',bsfp(1),bsfp(2),bsfp(6),bsfp(7)))
xticklabels(ax3,{})
ax3 = gca; 
ax3.FontSize = 16; 

ax4 = nexttile;
plot(t,SSTforced,'Color',[1 0.6 0.2],LineWidth = 5)
xlabel('Time (hrs)') 
ylabel('SST (^oC)') 
ylim([min(SSTforced)-0.5, max(SSTforced)+0.5])
xlim([tspan(1),tspan(2)])
% title('SST forcing profile',...
%     sprintf('SST = [%.1f,%.1f] %cC, rate = %.2f, peak = %d hrs.',Ssfp(1)-273.15,Ssfp(2)-273.15,char(176),Ssfp(5),Ssfp(6)))
ax4 = gca; 
ax4.FontSize = 16; 

til.TileSpacing = 'compact';

%% plot stacked
% Irma wind data
load('IrmaWindData.mat')

% colours
blu = [0 128 255]./255;
lblu = [0 255 255]./255;
gre = [0 255 51]./255;
red = [255 51 51]./255;
lred = [255 153 153]./255;
br = [160 82 45]./255;
dg = [77 149 66]./225;
gra = [160 160 160]./225;
pur = [255 0 255]./255;
ora = [255 128 0]./255;

figure(3); clf; hold on;
til = tiledlayout(2,1);

% tangential wind and radius
ax1 = nexttile;
yyaxis left; hold on;
plot(time,IrmaWind,'-.','Color',lred,LineWidth = 6) % 1:11, 17:27 
plot(t(1:24:end),windy(1:24:end),'-','Color',red,LineWidth = 6)
ylabel({'{\it v_{b2}} (ms^{-1})'})
ylim([min(windy)-5, max(windy)+5])
ax1.YColor = red;

yyaxis right
plot(t(1:24:end),rad(1:24:end),':','Color',blu,LineWidth = 6)
ylabel({'{\it r_{b2}}'; '(\times 10^{4}m)'}) 
ylim([min(rad)-0.1, max(rad)+0.1])
ax1.YColor = blu;

xlim([tspan(1),tspan(2)]) 
xticklabels(ax1,{})
ax1 = gca; 
ax1.FontName = 'Calibri';
ax1.FontSize = 22; 
box on;

% beta and SST
ax2 = nexttile;
yyaxis left
plot(t,SSTforced,'Color',dg,LineWidth = 6)
ylabel('SST (^oC)')
ylim([min(SSTforced)-0.1, max(SSTforced)+0.1])
ax2.YColor = dg;

yyaxis right
plot(t,betaforced,':','Color',pur,LineWidth = 6)
ylabel('\beta') 
ylim([min(betaforced)-0.05, max(betaforced)+0.05])
ax2.YColor = pur;

xlabel('Time (hrs)')
xlim([tspan(1),tspan(2)]) 
ax2 = gca; 
ax2.FontName = 'Calibri';
ax2.FontSize = 22;
box on;

til.TileSpacing = 'compact';

%% plot over SN boundary
load('SSTbetaSNCont.mat','SST_sn31','beta_sn31');

figure(4); clf;  hold on;
plot(beta_sn31,SST_sn31,linewidth = 3);
plot(betaforced,SSTforced,linewidth = 3)





