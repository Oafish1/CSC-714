% General Purpose
% Choose p (n,h)
p = 5;
n = 2^p - 1;
h = 1/(n+1);

% Initialize U
U = zeros(n*(n+2), 1);

% Initialize residuals
eff = zeros(n*(n+2), 1);
% Dirichlet BC
for i = 1:n+2
    eff(i) = f((i-1) * h);
end

% Apply multigrid method
smooth = 40;
U = multigrid(U, eff, smooth, 2^2);
U = reshape(U, [n+2, n]);
% Add Dirichlet boundaries for viewing
U_in = U;
U = zeros(n+2);
for i = 1:(n+2)
    U(i, 1) = f(h * (i-1));
end
for i = 1:n
    U(:, 1+i) = U_in(:, i);
end

% Show preview
sa = figure('Name','Multigrid Approximation','NumberTitle','off');
X = 0:h:1;
Y = 0:h:1;
surf(X, Y, U);
title('Multigrid: ' + string(h));


% Compare to Regular Method
% Choose p (n,h)
p = 5;
n = 2^p - 1;
h = 1/(n+1);

% Initialize U
U = zeros(n*(n+2), 1);

% Initialize Residuals
eff = zeros(n*(n+2), 1);
% Dirichlet BC
for i = 1:n+2
    eff(i) = f((i-1) * h);
end

% Calculate Fine Answer
p_star = p;
n_star = 2^p_star - 1;
h_star = 1/(n_star+1);

eff_star = zeros(n_star*(n_star+2), 1);
for i = 1:n_star+2
    eff_star(i) = f((i-1) * h_star);
end

U_star = zeros(n_star*(n_star+2), 1);
U_star = multigrid(U_star, eff_star, 70, 2^2);

for sh = 1:(p_star - p)
    U_star = shrink(U_star);
end

'Fine approximation created'

% Set testing vars
smooth = 51;
step = 1;

% Set vars
updates_per_multi_smooth = 0;
for power = 2:p %approximate
    updates_per_multi_smooth = updates_per_multi_smooth + n*(n+2)/(4^(p-power));
end

% Initialize placeholders
err_hist = zeros(2, ((smooth-1)/step) + 1);
U_gs = U;

for sm = 1:step:smooth
    % Apply multigrid method
    U_multi = multigrid(U, eff, sm, 2^2);
    err_hist(1, ((sm-1)/step) + 1) = 2*sm*updates_per_multi_smooth;
    err_hist(2, ((sm-1)/step) + 1) = norm(U_multi - U_star, inf);
    
    iter = 3;
    for it = 1:(iter*step)
        U_gs = jacobi(U_gs, eff, 1);
    end
    err_hist(3, ((sm-1)/step) + 1) = (((sm-1)/step) + 1)*iter*step*(n*(n+2));
    err_hist(4, ((sm-1)/step) + 1) = norm(U_gs - U_star, inf);
end

% Show preview
sa = figure('Name','Error: Multigrid vs Jacobi','NumberTitle','off');
err_hist
plot(err_hist(1, :), err_hist(2, :))
hold on;
plot(err_hist(3, :), err_hist(4, :))
hold off;
title('Multigrid: ' + string(h));
xlabel('Grid Value Updates');
ylabel('Error');
%ylim([0, 1]);
legend('Multigrid', 'Jacobi');


% Main V-cycle multigrid implementation
function U = multigrid(U, eff, smooth, stop)
    w = 2/3;
    
    % Extract n
    U_size = size(U);
    n = sqrt(U_size(1)+1) - 1;
    
    % Extract h (Assumes Omega \in [0, 1]^2)
    h = 1/(n+1);
    
    % Apply simple iterative method 'smooth' times
    for i = 1:smooth
        U = jacobi(U, eff, w);
    end
    
    % Estimate ~A~e = ~-r (assumes that n+1 is a power of 2)
    r = residual(U, eff);
    r_down = shrink(r);
    
    err = zeros(size(r_down));
    if n > stop
        err = multigrid(err, r_down, smooth, stop); % No -r because C is already inverted
    end % Can also solve in an else statement here
    
    % Scale ~e
    err_scaled = interpolate(err);
    
    % Update U with ~e
    U = U + err_scaled; % Notice +, this is because we take -~\Delta
    
    % Apply simple iterative method 'smooth' times
    for i = 1:smooth
        U = jacobi(U, eff, w);
    end
end

function shrunk = shrink(orig)
    % Extract n
    orig_size = size(orig);
    n = sqrt(orig_size(1)+1) - 1;
    
    orig = reshape(orig, [n+2, n]);
    shrunk = orig(1:2:(n+2), 2:2:(n-1)); % Note that x boundaries are not stored in U
    shrunk = reshape(shrunk, 1, []).';
end

function scaled = interpolate(orig)
    % Extract n
    orig_size = size(orig);
    n = sqrt(orig_size(1)+1) - 1;
    
    orig = reshape(orig, [n+2, n]);
    
    scaled = zeros(n*2 + 3, n*2 + 1);
    size_orig = size(orig);
    
    % Set base
    for i = 1:size_orig(1)
        for j = 1:size_orig(2)
            scaled(2*(i-1) + 1, 2*j) = orig(i, j);
        end
    end
    
    % Interpolate (could be replaced with clone to save time)
    % y
    for i = 1:(size_orig(1)-1)
        scaled(2*i, :) = (scaled(2*i - 1, :) + scaled(2*i + 1, :))/2;
    end
    % x boundaries (clone)
    scaled(:, 1) = scaled(:, 2);
    scaled(:, n*2+1) = scaled(:, n*2);
    % x
    for j = 2:size_orig(2)
        scaled(:, 2*(j-1) + 1) = (scaled(:, 2*(j-1)) + scaled(:, 2*j))/2;
    end
    
    scaled = reshape(scaled, 1, []).';
end

function r = residual(U, eff)
    % Extract n
    U_size = size(U);
    n = sqrt(U_size(1)+1) - 1;
    
    % Get A
    C = CMatrix(n);
    
    % Calculate R
    r = eff - C * U;
end

% Single iteration of the Jacobi algorithm
function U = jacobi(U, r, w)
    % Extract n
    U_size = size(U);
    n = sqrt(U_size(1)+1) - 1;
    
    % Extract h (Assumes Omega \in [0, 1]^2)
    h = 1/(n+1);
    
    % Iterate (Assumes 2D Poisson)
    U = (1-w) * U + w * (1/4) * ((4*eye(n*(n+2)) - CMatrix(n)) * U + r);
end

function C = CMatrix(n)
    % A
    A = zeros(n, n);
    A(1, 2) = -1;
    A(n, n-1) = -1;
    for i = 2:(n-1)
        A(i, i+1) = -1;
        A(i, i-1) = -1;
    end
    
    % B
    B = zeros(n+2, n+2);
    % First and last rows
    B(1, 1) = 4;
    B(1, 2) = -2; % To account for Neumann BC
    B(n+2, n+1) = -2; % To account for Neumann BC
    B(n+2, n+2) = 4;
    for i = 2:n+1
        B(i, i+1) = -1;
        B(i, i) = 4;
        B(i, i-1) = -1;
    end
    
    % C, formally known as A
    C = kron(A, eye(n+2)) + kron(eye(n), B);
end

% Dirichlet Boundary function
function out = f(y)
    out = cos(2 * pi * y);
    %out = sign(cos(2 * pi * y));
end