%% sea surface entropy under ambient region
function s_ocean_a = S_s_oa(s_i,p)

    s_oa0 = p(40);

    s_ocean_a = (S_s_oi(s_i,p)-s_oa0)/2;
    end
