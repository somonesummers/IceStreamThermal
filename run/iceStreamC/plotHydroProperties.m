clear all
clc

%%%%%     Plot Discharge     %%%%%
box on
%Coefficients defined by Perol et al. (2015)
n_m = 0.025, A_melt = 24e-25, K_2 = 0.05969, f = 0.4, alpha = 0.0015;
tau_max = [];
load data
tau_max = (3.7e3+76.3e3*(t/1e9).^2);
Q = (n_m^(1/4)*A_melt^(1/3)*tau_max./(K_2*f*alpha^(11/24))).^(1/12);
subplot(2,1,1)
plot(t/max(t),Q,'k','LineWidth',4)
xlim([0 1.0])
ylim([0.72 0.96])
ylabel('Discharge [m^3/s]')

%%%%%     Plot Film Thickness     %%%%%
channel_h = [];
for i = 1:length(t)
    [tau_E,epsilon_E] = effectiveStress(u(:,i),mu(:,i),m,n,dy,dz);
    b_dot = meltwaterFlux(m,nT,dz,T(:,i),T_m,[zeros(m*(nT-n),1);tau_E],...
            [zeros(m*(nT-n),1);epsilon_E],G_base,rho,L,u(:,i),tau_base(:,i)); %[m/s]
    h = (((12*1.787e-3)*(b_dot)/0.4)./gradient(gradient(tau_base(:,i),dy),dy)).^(1/3);
    channel_h = [channel_h,h(y==33200)];
end
subplot(2,1,2)
plot(t/max(t),channel_h*1e3,'k','LineWidth',4)

ylabel('Film thickness [mm]')
xlabel('Normalized time (t/T)')