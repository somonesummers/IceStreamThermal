clear all
clc

y = (0:1/100:1)';
u =@(y) (1-(1-y).^4);

%%%   Plot digitized velocities   %%%
load digitizedVelocities
% Ice Stream B (dnB)
% G = diag(u(dnB(:,1)/max(dnB(:,1))));
% E = (G'*G)\G'*(dnB(:,2)/max(dnB(:,2)));
subplot(2,4,5)
E2 = 1 + 0.2*smearedHeavi(0.8-y,0.8) + 0.75*smearedHeavi(0.38-y,0.4) + ...
     0.02*gaussmf(y,[0.04 0.32]);
plot(y,E2,'k','LineWidth',1.5)
% hold on
% plot(dnB(:,1)/max(dnB(:,1)),E)
xlabel('y/Y')
ylabel('E')
subplot(2,4,1)
plot(dnB(:,1)/max(dnB(:,1)),dnB(:,2)/max(dnB(:,2)),'.r','MarkerSize',20)
hold on
plot(y,E2.*u(y),'k','LineWidth',1.5)
% plot(dnB(:,1)/max(dnB(:,1)),E.*u(dnB(:,1)/max(dnB(:,1))))
axis([0 1 0 1])
title('Downstream Whillans')
ylabel('u/U')
% Ice Stream B (upB)
% G = diag(u(upB(:,1)/max(upB(:,1))-.05));
% E = (G'*G)\G'*(upB(:,2)/max(upB(:,2)));
subplot(2,4,6)
E2 = 1 + 0.1*smearedHeavi(0.8-y,0.8) + 0.44*smearedHeavi(0.5-y,0.5);
plot(y,E2,'k','LineWidth',1.5)
% hold on
% plot(upB(:,1)/max(upB(:,1))-.05,E)
xlabel('y/Y')
subplot(2,4,2)
plot(upB(:,1)/max(upB(:,1))-.05,upB(:,2)/max(upB(:,2)),'.r','MarkerSize',20)
hold on
plot(y,E2.*u(y),'k','LineWidth',1.5)
% plot(upB(:,1)/max(upB(:,1))-.05,E.*u(upB(:,1)/max(upB(:,1))-.05))
axis([0 1 0 1])
title('Upstream Whillans')
% Ice Stream E
% G = diag(u(ISE(:,1)/max(ISE(:,1))));
% E = (G'*G)\G'*(ISE(:,2)/max(ISE(:,2)));
subplot(2,4,7)
E2 = 1 + 0.4*smearedHeavi(0.17-y,0.17) + 0.9*smearedHeavi(0.35-y,0.35) + ...
     0.05*gaussmf(y,[0.1 0.3]) - 0.19*gaussmf(y,[0.15 0.5]);
plot(y,E2,'k','LineWidth',1.5)
% hold on
% plot(ISE(:,1)/max(ISE(:,1)),E)
xlabel('y/Y')
subplot(2,4,3)
plot(ISE(:,1)/max(ISE(:,1)),ISE(:,2)/max(ISE(:,2)),'.r','MarkerSize',20)
hold on
plot(y,E2.*u(y),'k','LineWidth',1.5)
% plot(ISE(:,1)/max(ISE(:,1)),E.*u(ISE(:,1)/max(ISE(:,1))))
axis([0 1 0 1])
title('MacAyeal')
% Ice Stream D
% G = diag(u(ISD_profileF(:,1)/max(ISD_profileF(:,1))));
% E = (G'*G)\G'*(ISD_profileF(:,2)/max(ISD_profileF(:,2)));
subplot(2,4,8)'
E2 = 0.99 + 0.35*smearedHeavi(0.45-y,0.45) - 0.02*gaussmf(y,[0.05 0.4]) + ...
     0.02*gaussmf(y,[0.05 0.3]);
plot(y,E2,'k','LineWidth',1.5)
% hold on
% plot(ISD_profileF(:,1)/max(ISD_profileF(:,1)),E)
xlabel('y/Y')
subplot(2,4,4)
plot(ISD_profileF(:,1)/max(ISD_profileF(:,1)),...
     ISD_profileF(:,2)/max(ISD_profileF(:,2)),'.r','MarkerSize',20)
hold on
plot(y,E2.*u(y),'k','LineWidth',1.5)
% plot(ISD_profileF(:,1)/max(ISD_profileF(:,1)),...
%      E.*u(ISD_profileF(:,1)/max(ISD_profileF(:,1))))
axis([0 1 0 1])
title('Bindschadler')
setFontSize(16)

% plot(y/1e3,u(y)*3.1557e7)
% axis([0 65 0 600])
% xlabel('Distance from margin [km]')
% ylabel('Velocity [m/yr]')
% hold on
% plot(B(:,1)/1e3,B(:,2))