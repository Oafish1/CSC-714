% Preset variables
N = 100;
Nt = N;
len = 10; % Extend observed time

% Initialize our matrix
u = zeros(N+1, N+1, Nt+1);
c = .2;
% Neumann boundary condition
ut0 = zeros(N+1);
for x = 2:N
    for y = 2:N
        ut0(x,y) = u(x,y,2) - f( (x-1) * (1/N) ) * ...
            f( (y-1) * (1/N) ) * (2/Nt);
    end
end
u(:,:,2) = 2 * u(:,:,1) - ut0 + (1/Nt)^2 * c^2 * Deltau( u(:,:,1) );

% Iterate the matrix
for t = 2:len*Nt
    u(:,:,t+1) = 2 * u(:,:,t) - u(:,:,t-1) + ...
        (1/Nt)^2 * c^2 * Deltau( u(:,:,t) );
end


speed = 2;
% Plot over time
for t = 1:speed:len*Nt+1
    surf( 0:(1/N):1, 0:(1/N):1, u(:,:,t), 'EdgeColor', 'none' );
    axis([0 1 0 1 -1 1]);
    drawnow
end



% Prerequisite functions
function r = Deltau(u)
    u_size = size(u);
    N = u_size(1) - 1;
    
    ddx = zeros(size(u));
    ddy = zeros(size(u));
    
%     % Boundary data points (Pseudo)
%     ddx(1,:) = N^2 * ( ( 2 * u(1,:) - u(2,:) ) - 2 * u(1,:) + ...
%         u(2,:) );
%     ddx(N+1,:) = N^2 * ( ( 2 * u(N+1,:) - u(N,:) ) - 2 * u(N+1,:) + ...
%         u(N,:) );
% 
%     ddy(:,1) = N^2 * ( ( 2 * u(:,1) - u(:,2) ) - 2 * u(:,1) + ...
%         u(:,2) );
%     ddy(:,N+1) = N^2 * ( ( 2 * u(:,N+1) - u(:,N) ) - 2 * u(:,N+1) + ...
%         u(:,N) );
%     
%     % Boundary data points (More accurate, but only works with
%     % smaller c)
%     ddx(1,:) = N^3 * ( 2 * u(1,:) - 5 * u(2,:) + ...
%         4 * u(3,:) - u(4,:) );
%     ddx(N+1,:) = N^3 * ( 2 * u(N+1,:) - 5 * u(N,:) + ...
%         4 * u(N-1,:) - u(N-2,:) );
% 
%     ddy(:,1) = N^3 * ( 2 * u(:,1) - 5 * u(:,2) + ...
%         4 * u(:,3) - u(:,4) );
%     ddy(:,N+1) = N^3 * ( 2 * u(:,N+1) - 5 * u(:,N) + ...
%         4 * u(:,N-1) - u(:,N-2) );

    % Inner data points
    for x = 2:N
        for y = 1:N
            ddx(x,y) = N^2 * ( u(x-1,y) - ...
                2 * u(x,y) + u(x+1,y) );
        end
    end
    for x = 1:N
        for y = 2:N
            ddy(x,y) = N^2 * ( u(x,y-1) - ...
                2 * u(x,y) + u(x,y+1) );
        end
    end

    % Calculate laplacian
    r = ddx + ddy;
end

function r = f(x)
    r = exp( -400 * ( x - .5 )^2 );
%     r = sin( 5 * pi * x ) * cos( 5 * pi * x )^2 + x^2 - .1 * exp(x);
end