function [T,k,c] = thermalRelax(m,n,nT,y,zT,dy,dz,T,k,k1,k2,c,c1,c2,...
                                tau_E,epsilon_E,rho,w,v,northT,southT,...
                                eastT,westT,T_m,T_atm,G_base,u,T_prev,dt,...
                                tau_base,smearT,maxT,omegaT,tol,...
                                timeDepFlag,fID,MTP)
%Perform under-relaxation on the thermal model

for j = 1:maxT
    T_old = T;
    T = thermalModel(m,n,nT,y,dy,dz,T,T_prev,dt,k,c,tau_E,...
                     epsilon_E,rho,w,v,northT,southT,...
                     eastT,westT,T_m,T_atm,G_base,u(1:m),...
                     tau_base,smearT,timeDepFlag,MTP);
    T = omegaT*T + (1-omegaT)*T_old;
    k = k1*exp(-k2*1e-3.*T);
    c = c1 + c2.*T;
    fprintf(fID,'\t\t[%d]Temperature Residual: %f \n',j,norm((T_old-T)./T_old));
    resT = norm((T_old-T)./T_old);
    if resT < tol
        break
    end
%     plotTemperatureField(m,nT,y,zT,T,T_m)
%     getframe;
end

end

