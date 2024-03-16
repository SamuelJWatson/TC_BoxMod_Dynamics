function [a,b] = bisect(a,b,Basin1,Basin2,tol)
%BISECT Use this to bisect between two basins
%   a - point in Basin1
%   b - point in Basin2
%   Basin1 - full set of points in Basin1
%   Basin2 - full set of points in Basin2

dist = norm(a-b);
while dist>tol
    c = median([a;b],1);
    test1 = c-Basin1;
    min1 = min(vecnorm(test1,2,2));
    test2 = c-Basin2;
    min2 = min(vecnorm(test2,2,2));
    if min1 < min2
        a=c;
    else
        b=c;
    end
    dist = norm(a-b);
end

end