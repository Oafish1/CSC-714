% To use, set the variables in the next line.  'sim_time' is to what t the
% simulation should run.  'simulate' defines whether or not there should be
% a live preview of the calculations.  The calculation method (Euler vs
% Godunov) can be set with the 'godunov' variable in solveIB().  Initial
% conditions can be set in f().  'ran' defines the y range of the
% simulation graphs.
n=100;sim_time=4;simulate=0;ran=[-3 3];
n_t=4*n;U=zeros(sim_time*n_t+1,n+1);
U(1,:)=f(0:1/n:1);
% Calculate Solution
fig=figure('Name','Simulation','NumberTitle','off');
for t=1:sim_time*n_t
    U(t+1,:)=solveIB(U(t,:),n_t);
    if simulate
        plot(0:1/n:1,U(t,:))
        title(strcat('t=',string(round(t/n_t,2))))
        xlabel('x')
        axis([0 1 ran]);
        drawnow
    end
end
close(fig)
% Plot with t as a Spatial Dimension
figure('Name','Summary','NumberTitle','off');
for s=1:sim_time
    subplot(1,sim_time,s)
    surf(0:1/n:1,0:1/n_t:sim_time,U)
    title(strcat('t=',int2str(s-1),'-',int2str(s)))
    xlabel('x')
    ylabel('t')
    axis([0 1 s-1 s ran]);
    colormap gray
    shading interp
end
% Plot x at Specified Intervals
figure('Name','Samples','NumberTitle','off');
samples=[0 .25 .5 1.5 3.25];
for s=samples
    subplot(1,size(samples,2),find(samples==s))
    t=s*n_t+1;plot(0:1/n:1,U(t,:))
    title(strcat('t=',string(sprintf('%.2f',t/n_t))))
    xlabel('x')
    axis([0 1 ran]);
    drawnow
end
% Calculate Error for Each h
figure('Name','Error','NumberTitle','off');
fine_n=1/.0003125;h=[.1 .05 .01 .005 .0025 .00125 .000625];err=[]; 
fine_n_t=4*fine_n;U_fine=zeros(fine_n_t+1,fine_n+1);U_fine(1,:)=f(0:1/fine_n:1);
for t=1:fine_n_t
    U_fine(t+1,:)=solveIB(U_fine(t,:),fine_n_t);
end
for n=h.^-1
    n_t=4*n;U=zeros(n_t+1,n+1);U(1,:)=f(0:1/n:1);
    for t=1:n_t
        U(t+1,:)=solveIB(U(t,:),n_t);
    end
    err=[err max(max(abs(...
        U-U_fine(1:fine_n_t/n_t:fine_n_t+1,1:fine_n/n:fine_n+1))))];
end
semilogx(h,err)
xlabel('h')
ylabel('Max Norm of Error')

% Solve with Periodic Boundary Conditions
function U_new=solveIB(U,n_t)
    conservative=0;godunov=1;
    U_new=U;n=size(U,2)-1;
    for x=1:n+1
        if godunov
            U_new(x)=U(x)-(n/n_t)*(F(U(mod(x,n)+1),U(mod(x-1,n)+1))...
                -F(U(mod(x-1,n)+1),U(mod(x-2,n)+1)));
        elseif conservative
            % u_t+(1/2)u^2_x=0
            u_x=(1/2)*n*(U(mod(x-1,n)+1)^2-U(mod(x-2,n)+1)^2);
            U_new(x)=U(x)-(1/n_t)*u_x;
        else
            % u_t+uu_x=0
            u_x=n*(U(mod(x-1,n)+1)-U(mod(x-2,n)+1));
            U_new(x)=U(x)-(1/n_t)*U(x)*u_x;
        end
    end
end
% Godunov Flux
function r=F(U_2,U_1)
    if U_1<=U_2
        if 0<=U_1
            r=U_1^2/2;
        elseif U_2<=0
            r=U_2^2/2;
        else
            r=0;
        end
    else
        if 0<=(U_2+U_1)/2
            r=U_1^2/2;
        else
            r=U_2^2/2;
        end
    end
end
% Initial Conditions
function r=f(x)
    r=3/2+sin(2*pi*x);
    %r=.5+sin(2*pi*x);
    %r=exp(-400*(x-.5).^2);
end