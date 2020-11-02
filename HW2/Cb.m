% Plot imaginary numbers by sample
X = -6:.01:2;
Y = -4:.01:4;

X_size = size(X);
Y_size = size(Y);

val = zeros(Y_size(2), X_size(2));

for x_in = 1:X_size(2)
    for y_in = 1:Y_size(2)
        x = X(x_in);
        y = Y(y_in);
        
        z = x + y * sqrt(-1);
        
        xi1 = abs( 2 + z + sqrt(z^2 + 4*z) ) / 2;
        xi2 = abs( 2 + z - sqrt(z^2 + 4*z) ) / 2;
        
        % Eliminate non-simple 1 roots
        if xi1 == xi2
            if xi1 == 1
                [x,y,xi1,xi2]
                xi2 = inf;
            end
        end
        
        % Shade if valid
        if xi1 <= 1
            if xi2 <=2
                val(y_in,x_in) = 1;
            end
        end
    end
end

surf(X,Y,val)
shading interp 
view(2);
