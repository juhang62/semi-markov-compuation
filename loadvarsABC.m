%load parameters and run computation
clear global
global kabs kbc0s kbas ksps delta1s kbt ksfs Lstep maxdis
Lstep=8; 
ksps=Lstep*[1 1]*0.34; %unit length = step size
maxdis=6; %6 or 3
kbt=1.38e-2*295;
F=-15:25;

%%%%kinetc parameters identical motors
kabs=5000*ones(1,2);
kbc0s=1000*ones(1,2);
kbas=0.2*kabs; ksfs=kabs./(0.0077*kabs-1);
delta1s=1*ones(1,2);
kdets=[5 5];
delta4=1; %detachment delta or 0.3

% %%kin1 and kin2
% kabs=[3000 3000]; kbc0s=[2753 1469];
% kbas=0.2*kabs; ksfs=[99 66.8];
% delta1s=[3 2];  %1 0.5
% kdets=[1 1];
% delta4=1; %detachment delta

% %%kin1 and kin1
% kabs=[3000 3000]; kbc0s=[2753 2753];
% kbas=0.2*kabs; ksfs=[99 99];
% delta1s=[3 3]; 
% kdets=[1 1];
% delta4=1; %detachment delta

%%%%%%%%%%%%
% deltas=0.62;
% katts=20;


% [Vsin1,Dsin1]=kinABCfor(F,1);
% [Vsin2,Dsin2]=kinABCfor(F,2);
% katt=50; 
% [V,D,stdistrcon,ER]=semimarkKineABCnew(F,katt,delta4,kdets); %katt,delta4,kdet1
[~,~,~,ER5]=semimarkKineABC(F,katt,delta4,[5 5]); %katt,delta4,kdet1
[~,~,~,ER10]=semimarkKineABC(F,katt,delta4,[20 20]);
[~,~,~,ER15]=semimarkKineABC(F,katt,delta4,[35 35]);
% figure
% plot(F,Vsin1,F,V)
% grid on
% title(['katt=' num2str(katt) '\delta=' num2str(delta4) 'k_{det}=' num2str(kdets(1))])


% %%%%%%%plot fcomp surf
% %deltas=0.1:0.1:1;
% deltas=linspace(0.1,2,10);
% %deltas=0.1
% %kdets=linspace(1,40,11);
% katts=linspace(50,200,11);
% %katts=5
% %katts=200;
% nkat=length(katts);
% ndel=length(deltas);
% %ndel=length(kdets);
% fcomp=zeros(ndel,nkat);
% [Vsin,Dsin]=kinABCfor(F,1);
% ind=find(Vsin<20,1);
% for ikat=1:nkat
%     for idel=1:ndel
%         delta4=deltas(idel);
%         %delta4=1;
%         %kdet=kdets(idel);
%         kdets=[5 5];
%         katt=katts(ikat);        
%         [V, D]=semimarkKineABCnew(F,katt,delta4,kdets);
%         Vdif=V-Vsin;
% %         [val, ind]=min(abs(Vdif(10:end)));
% %         ind=ind+9;
% %         fcomp(idel,ikat)=interp1(Vdif(ind-1:ind+1),F(ind-1:ind+1),0);
%         
%         fcomp(idel,ikat)=interp1(Vdif(16:ind),F(16:ind),0);
% %         if katts==125
% %             keyboard
% %         end
%     end
% end
% figure
% h=surf(katts,deltas,fcomp);
% xlabel('k_a_t_t')
% ylabel('\delta_4')
% zlabel('f_c_o_m_p')
% %title(['k_d_e_t=',num2str(kdets(1))])
% %zlim([2.6,3.2])
% view(45,45)