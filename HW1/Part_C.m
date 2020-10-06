% Subsection c
% Initialize vars
h_list = flip([.04, .1, .2, .5]);
err_hist = zeros(length(h_list), 1);
it_hist = zeros(length(h_list), 1);
max_iter = 100000;

% Make fine estimate
h_fine = .02; %must divide all in h_list
n_fine = (1/h_fine)-1;
[U_fine, ~, ~] = converge(h_fine, max_iter, zeros(n_fine+2), 10^(-13));

sa = figure('Name','Sample Approximations','NumberTitle','off');
% Test for multiple h
for h_i = 1:length(h_list)
    % Get the approximation
    h = h_list(h_i);
    [U, err, it] = converge(h, max_iter,...
        U_fine(1:int16(h/h_fine):int16(1/h_fine)+1,...
               1:int16(h/h_fine):int16(1/h_fine)+1),...
               10^(-9));
    err_hist(h_i) = err;
    it_hist(h_i) = it;
    
    % Show preview (fun)
    figure(sa);
    subplot(1, length(h_list), h_i);
    X = 0:h:1;
    Y = 0:h:1;
    surf(X, Y, U);
    title(h);
end

% % Show error
% fig = figure;
% X = 0:h:1;
% Y = 0:h:1;
% surf(X, Y, U - U_fine(1:int16(h/h_fine):int16(1/h_fine)+1,...
%                1:int16(h/h_fine):int16(1/h_fine)+1));

% Plot error
er = figure('Name','Error vs h','NumberTitle','off');
loglog(err_hist, h_list);
xlabel('Maximum Norm of Error');
ylabel('h');

% Plot iterations
it = figure('Name','h vs Iterations to Converge','NumberTitle','off');
loglog(h_list, it_hist)
xlabel('h');
ylabel('Iterations to Converge');


% Subsection f
% For problem E, only one member is used in h_list, h_fine.
% This was the only way I could find to avoid the differences discussed
% in C.g.  This is not necessary for C.f, only E.
% For example:
% h_list = flip([.01]);
% h_fine = .01;

% Initialize vars
h_list = flip([.01, .02, .05, .1]);
eps_list = [.5, .4, .3 .2, .1, .05, .02];
it_hist = zeros(length(h_list), length(eps_list));
it_hist_ex = zeros(length(h_list), length(eps_list));
max_iter = 20000;

% Make fine estimate
h_fine = .005; %must divide all in h_list
n_fine = (1/h_fine)-1;

% Initialize U_fine
U_fine = zeros(n_fine+2, n_fine+2);
% Dirichlet BC
for i = 1:n_fine+2
    U_fine(i, 1) = f((i-1) * h_fine);
end

% Create U_fine
% for i = 1:50000
%     U_fine = gaussSeidel(n_fine, U_fine);
% end
U_fine = converge(h_fine, max_iter, zeros(n_fine+2), 10^(-20));

% Test for multiple h
for h_i = 1:length(h_list)
    h = h_list(h_i);
    
    % Define n
    n = (1/h) - 1;

    % Initialize U
    it = 0;
    U = zeros(n+2, n+2);
    % Dirichlet BC
    for i = 1:n+2
        U(i, 1) = f((i-1) * h);
    end
    
    % Define U_fine with matching h
    U_fine_h = U_fine(1:int16(h/h_fine):int16(1/h_fine)+1, 1:int16(h/h_fine):int16(1/h_fine)+1);
    
    % Test for each epsilon (assumes list is sorted greatest to least)
    for eps_i = 1:length(eps_list)
        epsilon = eps_list(eps_i);
        
        % Begin iteration
        err = inf;
        while err > epsilon
            if it >= max_iter
                break
            end
            
            U = gaussSeidel(n, U);
            it = it + 1;
            
            err = norm(U - U_fine_h, inf);
        end
        
        it_hist(h_i, eps_i) = it;
        it_hist_ex(h_i, eps_i) = log(eps_list(eps_i))/log(...
            (1 - sin(pi/(2*(n+1)))^2 - sin(pi/(2*(n+3)))^2)^2 );
    end
end

% Only consider convergent cases (sets others to 0)
Z = (it_hist_ex - it_hist) .* (abs(it_hist) < max_iter);
% Raw iteration error
it = figure('Name','Iteration Prediction Error','NumberTitle','off');
X = eps_list;
Y = h_list;
if length(h_list) > 1
    surf(X, Y, Z);
    xlabel('epsilon');
    ylabel('h');
    zlabel('Prediction error');
else
    plot(X, Z);
    xlabel('epsilon');
    ylabel('Prediction error');
    title('IPE: h = ' + string(h_list(1)));
end

Z = (it_hist) .* (abs(it_hist) < max_iter) ./ it_hist_ex;
Z
% Fraction of expected iterations required
it = figure('Name','Iteration Prediction Fraction','NumberTitle','off');
X = eps_list;
Y = h_list;
if length(h_list) > 1
    surf(X, Y, Z);
    xlabel('epsilon');
    ylabel('h');
    zlabel('Fraction of expected iterations required');
else
    plot(X, Z);
    xlabel('epsilon');
    ylabel('Fraction of expected iterations required');
    title('IPF: h = ' + string(h_list(1)));
end


% Run Gauss-Seidel to convergence or max iteration
function [U, err, it] = converge(h, iter, fine, min_change_goal)
    % Define n
    n = (1/h) - 1;

    % Iterative approach
    U = zeros(n+2, n+2);
    % Dirichlet BC
    for i = 1:n+2
        U(i, 1) = f((i-1) * h);
    end

    % Begin iteration loop
    it = 0;
    change = inf;
    % Keep iterating as long as the error is sufficiently fluctuating
    while change > min_change_goal
        U_old = U;
        % Iterate by Gauss-Seidel
        U = gaussSeidel(n, U);
        
        % Break if over max iterations
        it = it + 1;
        if it >= iter
            break;
        end
        
        % Calculate error (for some norm)
        E = U - U_old;
        change = norm(E, inf);
    end
    % Calculate output error (for max norm, not particularly stable with the sign f(y))
    E = U - fine;
    err = norm(E, inf);
end

% Single iteration of the Gauss-Seidel algorithm
function U = gaussSeidel(n, U)
    % Neumann BC
    for j = 2:n+1
        U(1, j) = (1/4)*(2*U(2, j) + U(1, j+1) + U(1, j-1));
        U(n+2, j) = (1/4)*(2*U(n+1, j) + U(n+2, j+1) + U(n+2, j-1));
    end

    % Main iteration loop
    for i = 2:n+1
        for j = 2:n+1
            U(i, j) = (1/4)*(U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
        end
    end
end

% Single iteration of the Jacobi algorithm
function U_new = jacobi(n, U)
    U_new = U;
    % Neumann BC
    for j = 2:n+1
        U_new(1, j) = (1/4)*(2*U(2, j) + U(1, j+1) + U(1, j-1));
        U_new(n+2, j) = (1/4)*(2*U(n+1, j) + U(n+2, j+1) + U(n+2, j-1));
    end

    % Main iteration loop
    for i = 2:n+1
        for j = 2:n+1
            U_new(i, j) = (1/4)*(U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1));
        end
    end
end
    
function out = f(y)
    out = cos(2 * pi * y);
    %out = sign(cos(2 * pi * y));
end