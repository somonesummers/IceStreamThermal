function [M,b] = applyRobin(M,b,BC,dir,dx)
%Applies Robin BCs to stiffness matrix and source vector
%a*mu*(du/dn) + b*u = c
%BC: [[node numbers], [a], [mu], [b], [c]]
%dir = N: -m, S: +m, E: -1, W: +1

sortrows(BC,1);

M(BC(:,1),BC(:,1)) = M(BC(:,1),BC(:,1)) + ...
                     M(BC(:,1),BC(:,1)+dir).*speye(length(BC(:,1))).*...
                     -spdiags(2*dx*BC(:,4)./(BC(:,2).*BC(:,3)),0,...
                     length(BC(:,1)),length(BC(:,1)));

b(BC(:,1)) = b(BC(:,1)) - spdiags(M(BC(:,1),BC(:,1)+dir),0).* ...
             (2*dx*BC(:,5)./(BC(:,2).*BC(:,3)));
         
M(BC(:,1),BC(:,1)+dir) = M(BC(:,1),BC(:,1)+dir) + ...
                         M(BC(:,1),BC(:,1)+dir).*speye(length(BC(:,1)));

end