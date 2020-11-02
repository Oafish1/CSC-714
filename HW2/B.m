% Try for many N until error below threshold
for N = 20:1000
    % Sample the Function
    h = 1/N;
    samples = zeros(1, N+1);

    for j = 1:N+1
        samples(j) = f( j * h );
    end

    % Get the uniform norm
    diff = fm(0:.0000001:1) - interp1( (0:N)*h, fm( (0:N)*h ), 0:.0000001:1);
    
    % Note: epsilon = 10^-7 means a max error of 2 * 10^-5
    
    if max(diff) <= 10^-2 - 2 * 10^-5 % Calculated from E
        break
    end
end
[M, I] = max(diff);
fprintf( 'N = %d; x = %d; norm = %f\n', N, I, M )


function r = f(x)
    r = exp( -400 * ( x - .5 )^2 );
end

function r = fm(x)
    r = exp( -400 * ( x - .5 ).^2 );
end