% Plot imaginary numbers by sample
X = -6:.01:2;
Y = -2:.01:2;

X_size = size(X);
Y_size = size(Y);

val = zeros(Y_size(2), X_size(2));

for x_in = 1:X_size(2)
    for y_in = 1:Y_size(2)
        x = X(x_in);
        y = Y(y_in);
        
        % Make our complex number
        z = x + y * sqrt(-1);
        
        % Calculate complex roots
        xi1 = ( 2 + z + sqrt(z^2 + 4*z) ) / 2;
        xi2 = ( 2 + z - sqrt(z^2 + 4*z) ) / 2;
        
        % Eliminate non-simple 1 roots
        if xi1 == xi2 && abs(xi1) == 1
            [x,y,xi1,imag(xi1),xi2,imag(xi2)]
            xi2 = inf;
        end
        
        % Shade if valid
        if abs(xi1) <= 1 && abs(xi2) <= 1
            val(y_in,x_in) = 1;
        end
    end
end

surf(X,Y,val)
shading interp 
view(2);
