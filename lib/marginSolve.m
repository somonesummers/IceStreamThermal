function [u,T,mu,k,c,MTP] = marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,...
                                        u,T,mu,k,c,tau_base,enforcedMargin,...
                                        north,south,east,west,northT,southT,...
                                        eastT,westT,iceRheol,marginSolveType,...
                                        rho,g,alpha,drivingRamp,thermalEnhancement,...
                                        k1,k2,c1,c2,T_m,T_atm,G_base,w,v,...
                                        omegaM,omegaT,tol,timeDepFlag)
%Solve system of equations for self-consistent margin location

T_prev = T;
if marginSolveType == 1
    sweep = 0;
    MTP = sweep;
    iiRange = floor(log(Y/2)/log(2)):-1:0;
    jjRange = 1e10;
elseif marginSolveType == 2
    sweep = enforcedMargin-1;
    MTP = sweep;
    iiRange = 0;
    jjRange = 0;
end
for ii = iiRange
for jj = 0:jjRange
    offset = sweep + jj*2^ii; %offset: smearing region left position[m]
    eps = 2^ii; %eps: smearing region width[m]
    H_m = smearedHeavi(y-offset,eps); %H: smearedHeaviside function
    fprintf(fID,'Margin Position: %d \n',offset);
    fprintf(fID,'Margin Smearing Width: %d \n',eps);
    
    prevM = [u,mu];
    prevT = T;
    prevMTP = MTP;
    MTP = offset+eps;
    
    if iceRheol==1
        maxM = 1;
    elseif any(iceRheol==[2,3,4])
        maxM = 200;
    end
    maxT = 100;
    smearT = 4;
    for i = 1:maxM
        %Mechanical Model
        [u,tau_E,epsilon_E] = mechModel(m,n,y,dy,dz,Z,mu,tau_base,H_m,...
                                        north,south,east,west,rho,g,...
                                        alpha,drivingRamp);

        %Thermal Model
        if any(iceRheol==[3,4])
            [T,k,c] = thermalRelax(m,n,nT,y,zT,dy,dz,T,k,k1,k2,c,c1,c2,...
                                   tau_E,epsilon_E,rho,w,v,northT,southT,...
                                   eastT,westT,T_m,T_atm,G_base,u,T_prev,dt,...
                                   tau_base,smearT,maxT,omegaT,tol,...
                                   timeDepFlag,fID,MTP,thermalEnhancement);
        end
        
        %Rheological Model
        if any(iceRheol==[2,3,4,5])
            [mu,resM] = mechRelax(m,n,nT,z,T,mu,tau_E,rho,g,omegaM,iceRheol);
  	          fprintf(fID,'\t[%d]Viscosity Residual: %f \n',i,resM);
            if resM < tol
                break
            end
        end
        
        %Anomoly check for efficiency
        basal_dy = (mu(1:m)./(Z*rho*g*sin(alpha))).*gradient(u(1:m),dy);
        if (max(basal_dy) > 1e-1 && marginSolveType == 1)
            break
        end
    end
    
    %Anomoly check
    basal_dy = (mu(1:m)./(Z*rho*g*sin(alpha))).*gradient(u(1:m),dy);
    if (max(basal_dy(y >= offset)) > 1e-20 && marginSolveType == 1)
        sweep = offset;
        u = prevM(:,1);
        mu = prevM(:,2);
        T = prevT;
        k = k1*exp(-k2*1e-3.*T);
        c = c1 + c2.*T;
        MTP = prevMTP;
        fprintf(fID,'\n');
        break
    end

%     plotBasalStress(offset,eps,m,y,dz,Z,rho,g,alpha,tau_base,mu,u,marginPosition)
%     getframe;
end
end
fprintf(fID,'Final Margin Position: %.0f m \n',MTP);
[T,k,c] = thermalRelax(m,n,nT,y,zT,dy,dz,T,k,k1,k2,c,c1,c2,...
                       tau_E,epsilon_E,rho,w,v,northT,southT,...
                       eastT,westT,T_m,T_atm,G_base,u,T_prev,dt,...
                       tau_base,smearT,maxT,omegaT,tol,timeDepFlag,...
                       fID,MTP,thermalEnhancement);

end
