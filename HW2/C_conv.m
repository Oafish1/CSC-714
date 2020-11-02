% Preset variables
size = 5;
u_hist = zeros(1,size);
to_iter = [1000 500 250 200 100];
for it = 1:size
    N = to_iter(it);
    Nt = N;
    
    % Initialize our matrix
    u = zeros(N+1, N+1, Nt+1);
    c = .1;
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
    for t = 2:Nt
        u(:,:,t+1) = 2 * u(:,:,t) - u(:,:,t-1) + ...
            (1/Nt)^2 * c^2 * Deltau( u(:,:,t) );
    end
    
    if N == max(to_iter)
        u_fine = u;
    end
    
    u_hist(it) = max(max(max( ...
        u_fine( 1:( max(to_iter)/N ):max(to_iter)+1 , ...
        1:int16( max(to_iter)/N ):max(to_iter)+1 , ...
        1:int16( max(to_iter)/N ):max(to_iter)+1 ) - u )));
end

loglog(to_iter.^-1, u_hist)
title('Grid Spacing vs Error')
xlabel('h')
ylabel('Uniform Norm')




% Prerequisite functions
function r = Deltau(u)
    u_size = size(u);
    N = u_size(1) - 1;
    
    ddx = zeros(size(u));
    ddy = zeros(size(u));

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
end

function r = fm(x)
    r = exp( -400 * ( x - .5 ).^2 );
end