% Set Values
N = 24;

% Extend playtime len times
len = 5;

% Simulate
[U tic] = simulate(N,len);

% Plot over time
figure();
for t = 1:2:len*(N^2/6)+1
    %surf( tic, tic, U(:,:,t) );
    surf( tic, tic, U(:,:,t), 'EdgeColor', 'none' ); 
    axis([0 1 0 1 -1 1]);
    drawnow
end

% Error log plot
N_fine = 144;
N_list = [6 12 24 36 48];
% For less computation time
%N_fine = 48;
%N_list = [6 12 24];

% Calculate fine
w = waitbar(0,'Calculating fine solution...');
[U_fine tic] = simulate(N_fine,1);
waitbar(1/(1+length(N_list)),w,'Calculating solutions for smaller N...')

% Calculate error
err = 0;
for N_temp = N_list
    [U_temp tic] = simulate(N_temp,1);
    
    index = 1:(N_fine/N_temp):N_fine+1;
    index_t = 1:(N_fine^2/N_temp^2):N_fine^2/6+1;
    
    err = [err max( max( max( abs( U_fine(index,index,index_t)...
        - U_temp ))))];
    
    waitbar( ( 1+find(N_list==N_temp) ) / (1+length(N_list)),w )
end
close(w)

figure('Name','N vs Error','NumberTitle','off');
loglog([N_fine N_list], err);
xlabel('N');
ylabel('Maximum Norm of Error');


% Iteratively solve the heat equation
function [U tic] = simulate(N,len)
    % Calculated Values
    Nt = N^2 / 6;
    h = 1/N;
    ht = 1/Nt;

    % Chebyshev Grid Ticks
    tic_raw = cos( (0:N) * pi / N );
    tic = .5 + -( tic_raw / 2 );

    % Initialize grid
    U = zeros(N+1, N+1, Nt+1); % Delta t = Delta x,y

    % Solve for t = 1
    for it = 1:100
        % Iterate for t = 1
        U(:,:,2) = (1/4) * ( ...
            (25/12) * U(:,:,1) + ...
            (3) * U(:,:,3) - ...
            (4/3) * U(:,:,4) + ...
            (1/4) * U(:,:,5) + ...
            f(tic) .* f(tic)' * ht );

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
    for t = 2:len*Nt
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
function uxx = dxx(U, tic_raw, N)
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
    
    % Use differentiation matrix instead
    %uxx = dxxc(U,N);
end

% Second derivative using differentiation matrix
function uxx = dxxc(U,N)
    % Calculate using spectral FD scheme
    [ch ~] = cheb(N);
    uxx = U * ch' * ch';
    
    % Boundary conditions
    uxx(N+1,:) = 0;
    uxx(1,:) = 0;
    uxx(:,N+1) = 0;
    uxx(:,1) = 0;
end

% Compute differentiation matrix
% (From Trefethen, 'Spectral Methods in MATLAB')
function [D,x] = cheb(N)
    if N==0
        D=0;
        x=1;
        return
    end
    
    % Create the grid
    x = cos(pi*(0:N)/N)';
    X = repmat(x,1,N+1);
    
    % Calculate x_i - x_j
    dX = X-X';
    
    % Coefficients
    c = [2; ones(N-1,1); 2] .* (-1).^(0:N)';
    c = c * ( 1./c )';
    
    % Off-diagonal entries
    D  = c ./ ( dX + ( eye(N+1) ) );
    
    % Diagonals from identity 6.6
    D  = D - diag(sum(D'));
end

% Second derivative wrt y
function uyy = dyy(U, tic_raw, N)
    uyy = dxx( U',tic_raw,N )';
end


% Base function f
function r = f(x)
    r = exp( -400 * ( x - .4 ).^2 );
    %r = sin( 1*pi*x );
end