%% beta double sech forcing profile
% handed paramater vector with 1. minimum beta value, 2. maximum beta 
% value, 3. time of change between profiles, 4/8. decreasing (0) or 
% bincreasing (1) profile, 5/9. whether return (0)or ramp (1) profile, 
% 6/10. rate of change of forcing, 7/11. time of peak forcing. 
function beta = beta_sech(p, t)

betamin = p(1);
betamax = p(2);
change = p(3);

    if t <= change
        par = [ p(4);     %1 decreasing(0)/increasing(1) profile
                p(5);     %2 return(0)/ramp(1) profile
                p(6);     %3 rate
                p(7)];    %4 time of peak forcing
    elseif t > change
        par = [ p(8);     %1 decreasing(0)/increasing(1) profile
                p(9);     %2 return(0)/ramp(1) profile
                p(10);    %3 rate
                p(11)];   %4 time of peak forcing
    end
    
% i.e. (minimum beta) - (change in beta)x(forcing profile)
beta = betamin + (betamax - betamin)*forcing_sech(par,t);

end