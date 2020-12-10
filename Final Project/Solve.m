% Subroutine for solving Ax=f
function x=Solve(A,f)
    % Matlab Method Override
    %x=A\f;return;
    % Initialize
    max_iter=1000;
    x=zeros(size(f));x_old=ones(size(f));
    it=0;
    while norm(x-x_old,inf)>10^-12 && it<max_iter
        x_old=x;
        x=f-(A*x-diag(A).*x);
        x=x./diag(A);
        it=it+1;
    end
    if it==max_iter
        sprintf('Solve(): Max iterations exceeded')
    end
end