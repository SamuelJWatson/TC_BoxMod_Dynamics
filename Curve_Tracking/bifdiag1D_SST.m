function [ SN_lab,HB_lab,BP_lab ] = bifdiag1D_SST( bd,var,par )
%Plots 1D bifurcation diagram - COCO
%   bd - branch run ; var - variable index to plot ; par - continuation parameter
hold on

blu = [0 128 255]./255;
gre = [0 255 128]./255;
red = [255 51 51]./255;
br = [160 82 45]./255;
dg = [77 149 66]./225;

%extract columns
sst = coco_bd_col(bd,par);         % SST
x = coco_bd_col(bd,'x');         % entropy
ha = coco_bd_col(bd,'h_a');      % ambient region rel humidity
hb = coco_bd_col(bd,'h_refb');   % boundary layer rel humidity
bt = coco_bd_col(bd,'beta');     % exponent of radial decline

params = [sst; ha; hb; bt];
p = sst - 273.15.*ones(1,length(sst));

% convert entropy to tangential wind
windy = NaN(1,length(x(var,:)));
for i = 1:length(windy)
   windy(i) = v_b2(x(var,i),params(:,i)); 
end
xvar = windy;

%stability
stab = coco_bd_col(bd,'ep.test.USTAB');
changestab = [1];

for i = 2:length(stab)
    if stab(i) ~= stab(i-1)
        changestab = [changestab,i];
    end
end

for j = 1:length(changestab)-1
    ind_change = changestab(j);
    if stab(ind_change) == 0
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'color',blu,...
            'LineWidth',2);
    elseif stab(ind_change) == 1
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'color',gre,...
            'LineWidth',2);
    else
        plot(p(ind_change:changestab(j+1)),xvar(ind_change:changestab(j+1)),'color',red,...
            'LineWidth',2);
    end
end
       
ind_change = changestab(end);
if stab(ind_change) == 0
    plot(p(ind_change:end),xvar(ind_change:end),'color',blu,'LineWidth',2);
elseif stab(ind_change) == 1
    plot(p(ind_change:end),xvar(ind_change:end),'color',gre,'LineWidth',2);
else
    plot(p(ind_change:end),xvar(ind_change:end),'color',red,'LineWidth',2);
end

% Plot special points using indices
%extract labels
SN = coco_bd_idxs(bd,'SN');
HB = coco_bd_idxs(bd,'HB');
BP = coco_bd_idxs(bd,'BP');

SN_lab = coco_bd_labs(bd,'SN');
HB_lab = coco_bd_labs(bd,'HB');
BP_lab = coco_bd_labs(bd,'BP');

plot(p(SN),xvar(:,SN),'o','color',dg,'LineWidth',2,'MarkerSize',8);
plot(p(HB),xvar(:,HB),'d','color',br,'LineWidth',2,'MarkerSize',8);
plot(p(BP),xvar(:,BP),'ko','LineWidth',2,'MarkerSize',8);
end
