function [M,b] = applyDirichlet(M,b,BC)
%Applies Dirichlet BCs to stiffness matrix and source vector
%BC: [[node numbers], [nodal values]]

sortrows(BC,1);

M(BC(:,1),:) = sparse(length(BC(:,1)),length(M));
b = b - M(:,BC(:,1))*BC(:,2);
b(BC(:,1)) = BC(:,2);
M(:,BC(:,1)) = sparse(length(M),length(BC(:,1)));
M(BC(:,1),BC(:,1)) = speye(length(BC(:,1)),length(BC(:,1)));

end