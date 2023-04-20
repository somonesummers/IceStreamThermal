clear
m = 100;
n = 10;
T_range = linspace(253, 273,n)';
tau_range = logspace(2, 8, m);
mu_range = zeros(size(n,m));
mu_range_2 = zeros(size(n,m));

for j = 1:numel(T_range)
    for i = 1:numel(tau_range)
        mu_range(j,i) = iceRheology_suckale(T_range(j),tau_range(i),1,1,900,9.8,1);
    end
end

for j = 1:numel(T_range)
    for i = 1:numel(tau_range)
%         mu_range_glen(j,i) = iceRheology_suckale(T_range(j),tau_range(i),1,1,900,9.8,1);
         mu_range_2(j,i) = iceRheology_goldsby2000(T_range(j),tau_range(i),1,1,900,9.8,1);
    end
end


%%
figure(1)
clf
loglog(tau_range,mu_range,'-')
hold on
loglog(tau_range,mu_range_2,'--')
% loglog(tau_range,1e32*tau_range.^(-3),'k--')
% loglog(tau_range,1e23*tau_range.^(-1.4),'k--')
% loglog(tau_range,1e20*tau_range.^(-0.8),'k--')
% loglog(tau_range,1e17*tau_range.^(0),'k--')
% loglog(tau_range,1e23*tau_range.^(-2),'c-.')
xlabel('Tau stress Pa')
ylabel('mu viscosity')
% xlim([1e1,1e9])
% ylim([1e5,1e20])