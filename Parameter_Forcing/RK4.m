%% fourth order runge-kutta for system of ODEs 

function [t,x,windy,radius] = RK4(f,tspan,h,x0)
 
    t = tspan(1):h:tspan(2); % create timesteps 
    x = zeros(length(x0),length(t)); % create solution matrix
    x(:,1) = x0; % put inital values in
    windy = NaN(1,length(t)); % for wind conversion
    radius = NaN(1,length(t)); % for radius values
    
    % 4th order runge-kutta method
    for i = 1 : length(x)-1
        [K1,p] = f( t(i)      , x(:,i)          );  %get parameters as well
        [K2,~] = f( t(i) + h/2, x(:,i) + K1*h/2 );
        [K3,~] = f( t(i) + h/2, x(:,i) + K2*h/2 );
        [K4,~] = f( t(i) + h  , x(:,i) + K3*h   );
        x(:,i+1) = x(:,i) + (h/6)*( K1 + 2*K2 + 2*K3 + K4 );
        
        %convert to max tangential wind
        windy(i) = S_v_b2(x(1,i),p);
        
        %convert to radius
        radius(i) = S_r_b2(x(1,i),p);
    end
end