N_list = 20:10:100;

eig_max = [];
eig_min = [];
eig_min_4 = [];
eig_pred = [];
for N = N_list
    [ch x] = cheb(N);
    
    % Calculate D^2 and D^4
    ch = ch^2;
    ch_4 = 4 * ch(2:N,2:N)^2;
    ch = 2 * ch(2:N,2:N);
    
    % Store values
    eig_max = [eig_max max(eig(ch))];
    eig_min = [eig_min min(eig(ch))];
    eig_min_4 = [eig_min_4 min(eig(ch_4))];
    eig_pred = [eig_pred, -.096*N^4];
end

plot(N_list, [eig_max; eig_min_4; eig_min; eig_pred])
axis([min(N_list) max(N_list) -5*10^6 1*10^6])
title('Eigenvalues of 2D_N^2')
ylabel('Eigs')
xlabel('N')
legend('Max', 'Min 4D_N^4', 'Min', '-.096N^4')


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