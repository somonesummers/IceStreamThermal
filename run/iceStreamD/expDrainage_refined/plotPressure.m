clear all
clc

%%%%%     Plot Discharge     %%%%%
hold on
box on
n_m = 0.025;
A_melt = 24e-25;
K_2 = 0.05969;
f = 0.4;
alpha = 0.0015;
tau_max = [];
for i = 1:5
    cases = ['profileE';'profileD';'profileC';'profileB';'profileA'];
    cd(cases(i,:))
    load data
    tau_max = [max(-tau_base(y<=MTP)),tau_max];
    cd ..
end
dist_along = [0,15.0,27.7,40.6,53.4];
Q = (n_m^(1/4)*A_melt^(1/3)*tau_max./(K_2*f*alpha^(11/24))).^(1/12);
subplot(2,1,1)
plot(dist_along,Q,'k','LineWidth',4)
p = polyfit(dist_along(1:2),Q(1:2),1);
hold on
plot(dist_along,polyval(p,dist_along),'k--')
xlim([0 53.4])
ylim([0.81 0.9])
ylabel('Discharge [m^3/s]')
subplot(2,1,2)
%Using average b_dot
% plot(dist_along,[0.0022,0.0017,0.0012,8.09e-4,3.724e-4]*1e3,'k','LineWidth',4)
%Using actual b_dot
plot(dist_along,[9.49e-4,7.28e-4,4.71e-4,3.26e-4,1.65e-4]*1e3,'k','LineWidth',4)
axis([0 53.4 0 1.0])
ylabel('Film thickness [mm]')
xlabel('Distance along stream [km]')
hold off

%%%%%     Plot Transmissivity     %%%%%
% subplot(2,1,2)
% T_eff = [18.0e3,12.5e3,7.0e3,4.13e3,1.49e3];
% plot(dist_along,T_eff,'k','LineWidth',4)
% xlim([0 53.4])
% ylabel('Effective transmissivity [m^3/Pa s]')
% xlabel('Distance along stream [km]')

%%%%%     Plot film thickness     %%%%%
% subplot(2,1,2)
% plot(dist_along,(1./T_eff*12*1.787e-3).^(1/3),'k','LineWidth',4)
% xlim([0 53.4])
% ylabel('Effective transmissivity [m^3/Pa s]')
% xlabel('Distance along stream [km]')

%%%%%     Plot Meltwater Production     %%%%%
figure
meltwater = [];
for i = 2:6
    cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
    cd(cases(i,:))
    load data
    [tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz);
    b_dot = meltwaterFlux(m,nT,dz,T,T_m,[zeros(m*(nT-n),1);tau_E],...
            [zeros(m*(nT-n),1);epsilon_E],G_base,rho,L,u,tau_base); %[m/s]
    b_dot_int = dy*trapz(b_dot);
    meltwater = [meltwater,-b_dot_int];
    cd ..
end
hold on
plot(dist_along,fliplr(meltwater))
plot(dist_along,p(1)/1e3*ones(5,1))
ylabel('Meltwater production [m^2/s]')
xlabel('Distance along stream [km]')

%%%%%     Plot Actual Transmissivity     %%%%%
% figure
% meltwater = [];
% for i = 2:6
%     cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
%     cd(cases(i,:))
%     load data
%     [tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz);
%     b_dot = meltwaterFlux(m,nT,dz,T,T_m,[zeros(m*(nT-n),1);tau_E],...
%             [zeros(m*(nT-n),1);epsilon_E],G_base,rho,L,u,tau_base); %[m/s]
%     hold on
%     plot(y,(((12*1.787e-3)*(b_dot)/0.4)./gradient(gradient(tau_base,dy),dy)).^(1/3))
%     h = (((12*1.787e-3)*(b_dot)/0.4)./gradient(gradient(tau_base,dy),dy)).^(1/3);
%     h(y==MTP)
%     cd ..
% end