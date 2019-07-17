clear all
clc

%%%%%     Plot Discharge     %%%%%
box on
%Coefficients defined by Perol et al. (2015)
n_m = 0.025, A_melt = 24e-25, K_2 = 0.05969, f = 0.4, alpha = 0.0015;
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
xlim([0 53.4])
ylim([0.81 0.9])
ylabel('Discharge [m^3/s]')

%%%%%     Plot Film Thickness     %%%%%
subplot(2,1,2)
plot(dist_along,[9.49e-4,7.28e-4,4.71e-4,3.26e-4,1.65e-4]*1e3,'k','LineWidth',4)
axis([0 53.4 0 1.0])
ylabel('Film thickness [mm]')
xlabel('Distance along stream [km]')
hold off