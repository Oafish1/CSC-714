% The program takes longer than expected because of the dynamically
% created A().  This time is counted in delay and subtracted from the final
% total.
global cache
global delay
cache={};% Comment this line if not changing the initialization variables.
% Initialization
n_x=200;n_y=100;PML=3;J=100;c=1;
h=1/n_y;omega=1/h;k=omega/c;n_x_PML=n_x+2*PML+1;n_y_PML=n_y+2*PML+1;

% Define f
ax=[zeros(1,PML) 0:1/n_y:n_x/n_y zeros(1,PML)];
ay=[zeros(1,PML) 0:1/n_y:1 zeros(1,PML)];
f=sin(2*pi*ay)'*sin(2*pi*ax);
f=reshape(f,[],1);

% Rapid Decomposition Method
rTime=[];rResi=[];gm=1:6;
fprintf('GMRES\tTime\tResidual\n')
for gmres_iter=gm
    % Time
    delay=0;rapid=tic;
    [tu psi u]=RD(f,n_x,n_y,PML,J,gmres_iter,c);
    rapidTime=toc(rapid)-delay;
    rTime=[rTime rapidTime];
    
    % Precision
    resid_rapid=real(A(n_x_PML,n_y_PML,PML,k,h)*reshape(u,[],1)-f);
    resid_rapid=reshape(resid_rapid,n_y_PML,[]);
    resid_rapid=resid_rapid(PML+1:end-PML,PML+1:end-PML);
    resid_rapid=reshape(resid_rapid,[],1);
    rResi=[rResi max(abs(resid_rapid))];
    
    fprintf('%5d\t%.3f\t%.3e\n',...
        gmres_iter,...
        rapidTime,...
        max(abs(resid_rapid)))
end

% Plot
figure()
plot(gm,rTime)
xlabel('GMRES Iterations')
ylabel('Time')
figure()
semilogy(gm,rResi)
xlabel('GMRES Iterations')
ylabel('Residual')
drawnow

% Standard Method
delay=0;standard=tic;
us=Solve(A(n_x_PML,n_y_PML,PML,k,h),f);
standardTime=toc(standard)-delay;

% Comparison
fprintf('\nMethod\tTime\tResidual\n')
resid_standard=real(A(n_x_PML,n_y_PML,PML,k,h)*reshape(us,[],1)-f);
resid_standard=reshape(resid_standard,n_y_PML,[]);
resid_standard=resid_standard(PML+1:end-PML,PML+1:end-PML);
resid_standard=reshape(resid_standard,[],1);
fprintf('Stand\t%.3f\t%.3e\n',...
    standardTime,...
    max(abs(resid_standard)))

% Plot results
figure()
px=0-(PML/n_y):1/n_y:n_x/n_y+(PML/n_y);
py=0-(PML/n_y):1/n_y:1+(PML/n_y);

subplot(5,1,1)
surf(px,py,real(reshape(tu,n_y_PML,[])),'EdgeColor','none');
axis([0 n_x/n_y 0 1 -1 1])
title('$\tilde u$','Interpreter','latex','FontWeight','bold','FontSize',14)
view(2)

subplot(5,1,2)
surf(px,py,real(reshape(psi,n_y_PML,[])),'EdgeColor','none');
axis([0 n_x/n_y 0 1 -1 1])
title('$\psi$','Interpreter','latex','FontWeight','bold','FontSize',14)
view(2)

subplot(5,1,3)
surf(px,py,real(reshape(u,n_y_PML,[])),'EdgeColor','none');
axis([0 n_x/n_y 0 1 -1 1])
title('Rapid ($\tilde u+\psi$)','Interpreter','latex','FontWeight','bold','FontSize',14)
view(2)

subplot(5,1,4)
surf(px(PML+1:end-PML),py(PML+1:end-PML),...
    reshape(resid_rapid,n_y+1,[]),'EdgeColor','none');
axis([0 n_x/n_y 0 1 -1 1])
title('Rapid Residual ($Au-f$)','Interpreter','latex','FontWeight','bold','FontSize',14)
view(2)

subplot(5,1,5)
surf(px,py,real(reshape(us,n_y_PML,[])),'EdgeColor','none');
axis([0 n_x/n_y 0 1 -1 1])
title('Standard ($A \setminus f$)','Interpreter','latex','FontWeight','bold','FontSize',14)
view(2)

drawnow

function [tu psi u]=RD(f,n_x,n_y,PML,J,gmres_iter,c)
    % Initialize Workspace
    h=1/n_y;omega=1/h;k=omega/c;n_x_PML=n_x+2*PML+1;n_y_PML=n_y+2*PML+1;
    
    % Rapid Decomposition Method
    rapid=tic;
    % ag1
    % Define Boundaries
    b=PML+((0:J)/J)*n_x+1;
    cond=find((b-ceil(b)<.001).*([0 ones(1,J-1) 0]));
    b(cond)=b(cond)+.001;
    tb=b+[0 ones(1,J-1) 0];
    % ag1-

    % ag5
    % ag2
    % Calculate phi
    tu=bp(tb,f,n_y,PML,k,h);
    phi=f-A(n_x_PML,n_y_PML,PML,k,h)*reshape(tu,[],1);
    % ag2-

    % gmres/ag3
    psi_sub=gmres(phi,phi,...
        gmres_iter,0,@AP,b,tb,n_x,n_y,PML,k,h);
    % gmres/ag3-

    % ag4
    v=fp(b,psi_sub,n_y,PML,k,h);
    g=psi_sub-A(n_x_PML,n_y_PML,PML,k,h)*reshape(v,[],1); % Not used
    w=bp(tb,psi_sub,n_y,PML,k,h);
    psi=v+w;
    u=psi+tu;
    % ag4-
    % ag5-
end
% Backward pass
function tu=bp(tb,f,n_y,PML,k,h)
    n_y_PML=n_y+2*PML+1;dx=zeros(n_y_PML,1);f=reshape(f,n_y_PML,[]);
    tu=[];
    for j=size(tb,2)-1:-1:1
        fj=[zeros(n_y_PML,PML) f(:,ceil(tb(j)):floor(tb(j+1)))...
            zeros(n_y_PML,PML)];
        dd=[zeros(n_y_PML,PML)...
            zeros(n_y_PML,size(ceil(tb(j)):floor(tb(j+1))-1,2))...
            1/(2*h)*ones(n_y_PML,2) zeros(n_y_PML,PML-1)];
        right=reshape(fj+2*dd.*dx,[],1);
        tuj=Solve(A(size(fj,1),size(fj,2),PML,k,h),right);
        tuj=reshape(tuj,n_y_PML,[]);
        dx=tuj(:,2)-tuj(:,1);dx=dx/h; % asddf: PML effect?
        tuj=tuj(:,PML+1:end-PML);
        tu=[tuj tu];
    end
    tu=[zeros(n_y_PML,PML) tu zeros(n_y_PML,PML)];
end
% Forward pass
function tu=fp(b,f,n_y,PML,k,h)
    n_y_PML=n_y+2*PML+1;dx=zeros(n_y_PML,1);f=reshape(f,n_y_PML,[]);
    tu=[];
    for j=1:size(b,2)-1
        fj=[zeros(n_y_PML,PML) f(:,ceil(b(j)):floor(b(j+1)))...
            zeros(n_y_PML,PML)];
        dd=[zeros(n_y_PML,PML)...
            zeros(n_y_PML,size(ceil(b(j)):floor(b(j+1))-1,2))...
            1/(2*h)*ones(n_y_PML,2) zeros(n_y_PML,PML-1)];
        dd=fliplr(dd);
        right=reshape(fj+2*dd.*dx,[],1);
        tuj=Solve(A(size(fj,1),size(fj,2),PML,k,h),right);
        tuj=reshape(tuj,n_y_PML,[]);
        dx=tuj(:,end)-tuj(:,end-1);dx=dx/h; % asddf: PML effect?
        tuj=tuj(:,PML+1:end-PML);
        tu=[tu tuj];
    end
    tu=[zeros(n_y_PML,PML) tu zeros(n_y_PML,PML)];
end
% Bi-directional Pass (ag3)
function APf=AP(b,tb,f,n_x,n_y,PML,k,h)
    n_x_PML=n_x+2*PML+1;n_y_PML=n_y+2*PML+1;
    v=fp(b,f,n_y,PML,k,h);
    g=f-A(n_x_PML,n_y_PML,PML,k,h)*reshape(v,[],1);
    w=bp(tb,f,n_y,PML,k,h);
    he=g-A(n_x_PML,n_y_PML,PML,k,h)*reshape(w,[],1);
    APf=f-he;
end
% A^(j).  x_size is the size with PML for convenience.
function r=A(x_size,y_size,PML,k,h)
    global cache
    global delay
    create=tic;
    
    % Check if this A has already been created
    for a=cache
        if size(a{1},2)==x_size*y_size
            r=a{1};
            delay=delay+toc(create);
            return;
        end
    end
    
    % Laplacian w/ Homogenous Dirichlet
    Y=-4*eye(y_size)+diag(ones(1,y_size-1),1)+diag(ones(1,y_size-1),-1);
    X=diag(ones(1,x_size-1),1)+diag(ones(1,x_size-1),-1);
    L=kron(eye(x_size),Y)+kron(X,eye(y_size));
    L=L*h^2;

    % PML
    ax=[(1:-1/PML:1/PML).^2 zeros(1,x_size-2*PML) (1/PML:1/PML:1).^2];
    ax=(1+1i*ax).^-1;
    ay=[(1:-1/PML:1/PML).^2 zeros(1,y_size-2*PML) (1/PML:1/PML:1).^2];
    ay=(1+1i*ay).^-1;
    alpha=ay'*ax;
    L=L*diag(reshape(alpha,[],1)).^2;
    
    % Helmholtz
    r=L+k^2*eye(y_size*x_size);
    
    % Add to cache
    cache{end+1}=r;
    
    delay=delay+toc(create);
end
