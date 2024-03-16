%% SST Sech forcing profile 
% handed paramater vector with 1. minimum SST value, 2. maximum SST 
% value, 3. decreasing (0) or increasing (1) profile, 4. whether return (0)
% or ramp (1) profile, 5. rate of change of forcing, 6. time of peak
% forcing. 
function SST = SST_sech(p, t)

SSTmin = p(1);
SSTmax = p(2);

par = [ p(3);     %1 decreasing(0)/increasing(1) profile
        p(4);     %2 return(0)/ramp(1) profile
        p(5);     %3 rate
        p(6)];    %4 time of peak forcing

% i.e. (minimum SST) - (change in SST)x(forcing profile)
SST = SSTmin + (SSTmax - SSTmin)*forcing_sech(par,t);

end