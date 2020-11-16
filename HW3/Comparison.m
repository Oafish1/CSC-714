% Set Values
N_fine = 64;
N = 32;

% Calculate error
B_list = 5:5:40;

index = 1:(N_fine/N):N_fine+1;
index_f = 1 + (3/4) * ( N_fine^2 / 8 );
index_s = 1 + (3/4) * ( N^2 / 8 );

err_spec = [];
err_fini = [];

w = waitbar(0,'Calculating...');
for B = B_list
    % Calculate fine solutions
    % (Can also use interp, but that will favor the chosen method)
    [U_spec_fine ~] = simulate(N_fine,@dxx_spectral,@dyy_spectral,B,1);
    [U_fini_fine ~] = simulate(N_fine,@dxx_finite,@dyy_finite,B,0);
    
    % Calculate rough solutions
    [U_spec ~] = simulate(N,@dxx_spectral,@dyy_spectral,B,1);
    [U_fini ~] = simulate(N,@dxx_finite,@dyy_finite,B,0);
    
    % Record error
    err_spec = [err_spec log(max(max(max(abs(...
        U_spec_fine(index,index,index_f) - U_spec(:,:,index_s))))))];
    err_fini = [err_fini log(max(max(max(abs(...
        U_fini_fine(index,index,index_f) - U_fini(:,:,index_s))))))];
    
    % Update progress bar
    waitbar( find(B_list==B)/length(B_list),w )
end
close(w)

figure('Name','B vs Error at t=.75','NumberTitle','off');
plot(B_list, [err_spec; err_fini]);
legend('Spectral Error','Finite Error');
title('B vs Error at t=.75')
xlabel('B');
ylabel('log(Maximum Norm of Error)');


% Iteratively solve the heat equation
function [U tic] = simulate(N,dxx,dyy,B,grid)
    % Calculated Values
    Nt = N^2 / 8;
    h = 1/N;
    ht = 1/Nt;

    % Chebyshev Grid Ticks
    if grid
        tic_raw = cos( (0:N) * pi / N );
    else
        tic_raw = -1:(2*h):1;
    end
    tic = .5 + -( tic_raw / 2 );

    % Initialize grid
    U = zeros(N+1, N+1, Nt+1);

    % Solve for t = 1
    for it = 1:100
        % Iterate for t = 1
        U(:,:,2) = (1/4) * ( ...
            (25/12) * U(:,:,1) + ...
            (3) * U(:,:,3) - ...
            (4/3) * U(:,:,4) + ...
            (1/4) * U(:,:,5) + ...
            f(tic,B) .* f(tic,B)' * ht );

        % Solve for t = 2:4
        for t = 2:4
            % Second derivatives
            uxx = dxx( U(:,:,t),tic_raw,N );
            uyy = dyy( U(:,:,t),tic_raw,N );

            % Biharmonic
            uxxxx = dxx( uxx,tic_raw,N );
            uyyyy = dyy( uyy,tic_raw,N );
            uxxyy = dyy( uxx,tic_raw,N );
            
            % Solve for t
            U(:,:,t+1) = 2 * U(:,:,t) - U(:,:,t-1) + ...
                (ht)^2 * (uxx+uyy) + (ht^4/12) * (uxxxx+uyyyy+2*uxxyy);
        end
    end

    % Main loop
    for t = 2:Nt
        % Second derivatives
        uxx = dxx( U(:,:,t),tic_raw,N );
        uyy = dyy( U(:,:,t),tic_raw,N );

        % Biharmonic
        uxxxx = dxx( uxx,tic_raw,N );
        uyyyy = dyy( uyy,tic_raw,N );
        uxxyy = dyy( uxx,tic_raw,N );

        % Solve for t
        U(:,:,t+1) = 2 * U(:,:,t) - U(:,:,t-1) + ...
            (ht)^2 * (uxx+uyy) + (ht^4/12) * (uxxxx+uyyyy+2*uxxyy);
    end
end

% Second derivatives
% (Methodology from Trefethen, 'Spectral Methods in MATLAB')
function uxx = dxx_spectral(U, tic_raw, N)
    uxx = zeros(N+1);
    ii = 2:N;
    
    for x = 2:N
        % Choose our vector and make it periodic
        vec = U(x,:);
        periodic = [vec fliplr(vec(ii))];
        fourier = real(fft(periodic));
        
        % Make our vector V
        k_in = [0:N 1-N:-1];
        
        % Calculate (ik)^v
        coeff = ( 1i*k_in ).^1;
        % Set V_N = 0 -> W_N = 0 since v is odd
        coeff(N+1) = 0;
        % First derivative
        W1 = real(ifft( coeff .* fourier ));
        
        % Calculate (ik)^v
        coeff = ( 1i*k_in ).^2;
        % Second derivative
        W2 = real(ifft( coeff .* fourier ));
        
        % Second derivative by equation 8.7
        uxx(x,ii) = W2(ii) ./ ( 1 - tic_raw(ii).^2 ) - ...
            tic_raw(ii) .* W1(ii) ./ ( 1 - tic_raw(ii).^2 ).^(3/2);
    end
end

function uyy = dyy_spectral(U, tic_raw, N)
    uyy = dxx_spectral( U',tic_raw,N )';
end

% Finite second derivatives
function uxx = dxx_finite(U, ~, N)
    uxx = zeros(N+1);
    for y = 2:N
        uxx(:,y) = U(:,y+1) - 2*U(:,y) + U(:,y-1);
        uxx(:,y) = uxx(:,y) * N^2;
    end
    
    uxx(:,1) = 0;
    uxx(:,N+1) = 0;
    uxx(1,:) = 0;
    uxx(N+1,:) = 0;
end

function uyy = dyy_finite(U, tic_raw, N)
    uyy = dxx_finite( U',tic_raw,N )';
end

% Base function f
function r = f(x,B)
    r = sin( B*pi*x );
end