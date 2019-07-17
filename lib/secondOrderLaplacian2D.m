function [M] = secondOrderLaplacian2D(m,n,dx,dy,mu,a)
%Build stiffness matrix for Laplacian with non-uniform mu

%Reshape mu
mu_x = reshape(reshape(mu,m,n)',m*n,1);
mu_x = mean([[mu_x;mu_x(end-2*n+1:end-n)],[mu_x(n+1:2*n);mu_x]],2);
mu_y = mean([[mu;mu(end-2*m+1:end-m)],[mu(m+1:2*m);mu]],2);

%Build stiffness matrix
Aw = reshape(reshape(mu_x(1:end-n)/dx^2,n,m)',m*n,1);
Ae = reshape(reshape(mu_x(n+1:end)/dx^2,n,m)',m*n,1);
As = mu_y(1:end-m)/dy^2;
An = mu_y(m+1:end)/dy^2;
Ap = -(Ae+Aw+As+An);
M = spdiags([a.*An,a.*Ae,a.*Ap,a.*Aw,a.*As],[-m -1 0 1 m],m*n,m*n);

%Remove non-physical ghost point terms
M(m+1:m:m*(n-1)+1,m:m:m*(n-1)) = M(m+1:m:m*(n-1)+1,m:m:m*(n-1)) - ...
                                 M(m+1:m:m*(n-1)+1,m:m:m*(n-1)).*speye(n-1);
M(m:m:m*(n-1),m+1:m:m*(n-1)+1) = M(m:m:m*(n-1),m+1:m:m*(n-1)+1) - ...
                                 M(m:m:m*(n-1),m+1:m:m*(n-1)+1).*speye(n-1);

end