% Dimensions
N=400;Nt=1600;h=1/N;ht=1/Nt;PML=N/2;len=2;
% Initialization
U=zeros(1, 1+PML+N+1+PML+1);U_new=U;U_old=U;U_init=U;V=U;V_old=U;
% Initial Conditions
U_old(1+PML+1:1+PML+N+1) = 1*f((0:N)/N, .6);
% First Order Approximation wrt t
for x=2:size(U,2)-1
    V(x)=V_old(x)+(ht/h)*(U_old(x+1)-U_old(x-1));
end
for x=2:size(U,2)-1
    U(x)=U_old(x)+(ht/h)*(V(x+1)-V(x-1));
end
% Iterate
for t = 3:len*Nt+1
    V_new=zeros(1,1+PML+N+1+PML+1);
    for x=2:size(U,2)-1
        V_new(x)=V_old(x)+(ht/h)*(U(x+1)-U(x-1))-2*ht*sigma(x,N,PML)*V(x);
    end
    U_new=zeros(1,1+PML+N+1+PML+1);
    for x=2:size(U,2)-1
        U_new(x)=U_old(x)+(ht/h)*(V(x+1)-V(x-1))-2*ht*sigma(x,N,PML)*U(x);
    end
    V_old=V;V=V_new;U_old=U;U=U_new;
    % Plot
    dim=[-(PML+1)*h:h:0-h 0:h:1 1+h:h:1+(PML+1)*h];
    plot(dim,U);
    axis([min(dim) max(dim) -1 1]);
    line([0 0],[-1 1],'Color',[0 0 0])
    line([1 1],[-1 1],'Color',[0 0 0])
    drawnow
end
% Set Functions
function r=sigma(x,N,PML)
    scale=10;power=2;
    sig=0; if x<=1+PML, sig=(1+PML-x)/PML;
    elseif x>1+PML+N+1, sig=(x-(1+PML+N+2))/PML; end
    %if x<=1+PML||x>1+PML+N+1, sig=1; end
    sig=scale*sig^power;
    r=sig;
end
function r=f(x,c)
    r = exp(-400*(x-c).^2);
end