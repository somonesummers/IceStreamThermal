function [T] = thermalModel(m,n,nT,y,dy,dz,T,T_prev,dt,k,c,tau_E,epsilon_E,...
                            rho,w,v,northT,southT,eastT,westT,T_m,T_atm,...
                            G_base,u_base,tau_base,smearT,timeDepFlag,MTP)
%Solves thermal model

%%%%%     Property fields     %%%%%
G_base = G_base*ones(m,1);
tau_E = [zeros(m*(nT-n),1);tau_E];
epsilon_E = [zeros(m*(nT-n),1);epsilon_E];
H_t = smearedHeavi(T-T_m+smearT/2,smearT); %H_t: thermal Heaviside function

%%%%%     Different solution methods w/ and w/o till     %%%%%
if (n == nT) %w/o till
    %Build System
    M = secondOrderLaplacian2D(m,n,dy,dz,k,ones(m*n,1));
    b = -(1-H_t).*(2*tau_E.*epsilon_E); %b: shear heating[Pa/s]
%%%%%%%%%%     Note: change thermal heating factor here     %%%%%%%%%%
    b = 1*b;
    if (timeDepFlag == 1)
        M = speye(m*n,m*n) - (dt/rho)*secondOrderLaplacian2D(m,n,dy,dz,k,1./c);
        b = T_prev - (dt/(rho*c))*b;
    end
    
    %Apply BCs
    %North: stress-free, South: Robin, East,West: symmetry
    H_BC = heaviside(T(1:m)-T_m);
    H_BC(H_BC == 0.5) = 1;
    H_BC(1:find(H_BC==1,1)) = 1; %Enforce temperate to the left of temperate zone
%   H_BC = (y <= MTP); %Enforce temperate to the left of the MTP
    [M,b] = applyThermalBC(M,b,k,H_BC,G_base,T_atm,T_m,...
                           northT,southT,eastT,westT,m,n,nT,dy,dz);  
elseif (n ~= nT) %w/ till
    %Build System
    M = secondOrderLaplacian2D(m,nT,dy,dz,k,ones(m*nT,1));
    b = -(2*tau_E.*epsilon_E); %shear heating[Pa/s]
    b(m*(nT-n)+1:m+m*(nT-n)) = b(m*(nT-n)+1:m+m*(nT-n)) + ...
                               u_base.*tau_base; %frictional heating[Pa/s]
    b = (1-H_t).*b;
    if (timeDepFlag == 1)
        M = speye(m*nT,m*nT) - (dt/rho)*secondOrderLaplacian2D(m,nT,dy,dz,k,1./c);
        b = T_prev - (dt/(rho*c))*b;
    end
    
    %Apply BCs
    %North: stress-free, South: Neumann, East,West: symmetry
    [M,b] = applyThermalBC(M,b,k,1-H_t(1:m),G_base,T_atm,T_m,...
                           northT,southT,eastT,westT,m,n,nT,dy,dz);
end

%%%%%     Solve system     %%%%%
T = M\b;

end

% w = [zeros(m*(nT-n),1);reshape(repmat(w',m,1),m*n,1)];
% v = [zeros(m*(nT-n),1);reshape(repmat(v',m,1),m*n,1)];

% M = speye(m*nT,m*nT) - ...
%     (dt/(rho*c))*secondOrderLaplacian(m,nT,dy,dz,k) + ...
%     (dt)*secondOrderCentralAdvection(m,nT,dz,w,m) + ...
%     (dt)*secondOrderCentalAdvection(m,nT,dy,v,1);