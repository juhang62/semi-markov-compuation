%% 1-motor vs 2-motor k=500 vs 5000
set(0,'defaultAxesFontSize',14)
set(0, 'DefaultLineLineWidth', 2);
%subplot(1,3,1)
figure
plot(F,Vsin5h,'b',F,Vsin5t,'c--',F,V5h,'r',F,V5t,'m--')
xlim([min(F) max(F)])
xlabel('force (pN)')
ylabel('velocity (mm/s)')
%subplot(1,3,2)
figure
plot(F,Dsin5h,'b',F,Dsin5t,'c--',F,D5h,'r',F,D5t,'m--')
xlim([min(F) max(F)])
xlabel('force (pN)')
ylabel('diffusivity (mm^2/s)')
yticks(0:1000:4000)
legend('1-motor $k_{ab}$=500','1-motor $k_{ab}$=5000','2-motor $k_{ab}$=500','2-motor $k_{ab}$=5000','Interpreter','latex')
%subplot(1,3,3)
figure
plot(F,sqrt(Dsin5h)./Vsin5h,'b',F,sqrt(Dsin5t)./Vsin5t,'c--',F,sqrt(D5h)./V5h,'r',F,sqrt(D5t)./V5t,'m--')
xlim([min(F) max(F)])
xlabel('force (pN)')
ylabel('coef. of var.')

%% xs std bar plot
nkin=3;m=2;
xs=-maxdis:maxdis;
nxs=length(xs);
substd=stdistrcon(1:nkin^m*nxs,[11 16 21]); % F=10,0,10
eleM=ones(1,nkin^m);
contractM=blkdiag(eleM,eleM);
for i=3:nxs
    contractM=blkdiag(contractM,eleM);
end
xsstd=contractM*substd;
bar(xs,xsstd)
legend('F=-5pN','F=0','F=5pN')
ylim([0,0.4]);
xlabel('x_1-x_2')
ylabel('probablity')

%% kin1 kin2
plot(F,Vsin1,F,Vsin2,F,V)
xlabel('force (pN)')
ylabel('velocity (nm/s)')
legend('kin1','kin2','kin1-kin2')
figure
plot(F,Dsin1,F,Dsin2,F,D)
xlabel('force (pN)')
ylabel('diffusivity (nm^2/s)')
legend('kin1','kin2','kin1-kin2')

%% ER
semilogy(F,ER5,F,ER10,F,ER15)
xlabel('force (pN)')
ylabel('run length (nm)')
legend('k_{det}=5','k_{det}=20','k_{det}=35');

%% qss
plot(FF,V,FF,Vsemi,FFsim,Vsim,'o')
xlabel('force (pN)')
ylabel('velocity (nm/s)')
legend('qss1','qss2','MC')
plot(FF,D,FF,Dsemi,FFsim,Dsim,'o')
xlabel('force (pN)')
ylabel('diffusivity (nm^2/s)')
legend('qss1','qss2','MC')
