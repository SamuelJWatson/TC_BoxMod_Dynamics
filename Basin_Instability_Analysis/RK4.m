%% fourth order runge-kutta for system of ODEs 

function [t,x] = RK4(f,tspan,h,x0)
 
    t = tspan(1):h:tspan(2); % create timesteps 
    x = zeros(length(x0),length(t)); % create solution matrix
    x(:,1) = x0; % put inital values in
    
    % 4th order runge-kutta method
    for i = 1 : length(x)-1
        K1 = f( t(i)      , x(:,i)          );
        K2 = f( t(i) + h/2, x(:,i) + K1*h/2 );
        K3 = f( t(i) + h/2, x(:,i) + K2*h/2 );
        K4 = f( t(i) + h  , x(:,i) + K3*h   );
        x(:,i+1) = x(:,i) + (h/6)*( K1 + 2*K2 + 2*K3 + K4 );
    end
end