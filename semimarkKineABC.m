function [V, D, stdistrcon, ER]=semimarkKineABC(F,katt,delta4,kdets)
%semimarkKineABC  Given parameters, return asymptotic velocity, diffusivity
%   stationary distribution and expected run length using kinetic model ABC.

global  kabs kbc0s  kbas ksps delta1s kbt ksfs Lstep maxdis
%katt=20; 
%kdet1=1.11; %kdet2=7.4;
nkin=3; m=2; 
%%%% detach delta
delta3=0.3;%delta4=0.62;
%need to be continous at f=0 to have smooth curve
kdetfun=@(f,kdet1) kdet1*double(f>=0).*exp(f*delta3/kbt)+kdet1*double(f<0).*exp(-f*delta4/kbt); 
%kdetfun=@(f) 10;


%F=-15:15;
%F=-1:1;
%F=0;
nF=length(F);
V=zeros(nF,1); D=zeros(nF,1); ER=zeros(nF,1);

%%%%%%%
%xs=(-6:6)';
xs=(-maxdis:maxdis);
nxs=length(xs);
totst=nxs*nkin^m+6; %6 one-motor states; 7 if include 1 zero-motor state
stdistrcon=zeros(totst,nF);

%%%%%steping mapping
stptoP = [1 2 3];
stpfromP=[7 8 9];
stpZP=   [1 1 1]'/2;
stptoM = [1 4 7];
stpfromM=[3 6 9];
stpZM = [1 1 1]'/2;
detfrom=[1 2 3 1 4 7];
detto  =[1 2 3 4 5 6];

for iF=1:nF
ff=F(iF);
%%%%%faithful treatment of max displ but leads to bumpy curve
% xsmax=findxsmax(6,ff,ksps,maxdispl);
% xs=-xsmax:xsmax; %x1-x2!! depends on loading force
% xs=xs';

Q=zeros(totst);
Z1M=zeros(totst,1);Z2M=Z1M;

for i=1:nxs
    if i<nxs
        stpfromPi=stpfromP+(i-1)*nkin^m;
        stptoPi=stptoP+i*nkin^m;
        indtmp=sub2ind(size(Q),stpfromPi,stptoPi);
        indtmp=indtmp';
        Q(indtmp)=ksfs(1); 
        Z1M(stpfromPi)=Z1M(stpfromPi)+stpZP.*Q(indtmp);
        Z2M(stpfromPi)=Z2M(stpfromPi)+stpZP.^2.*Q(indtmp);
    end
    if i>1
        stpfromMi=stpfromM+(i-1)*nkin^m;
        stptoMi=stptoM+(i-2)*nkin^m;
        indtmp=sub2ind(size(Q),stpfromMi,stptoMi);
        indtmp=indtmp';
        Q(indtmp)=ksfs(2); 
        Z1M(stpfromMi)=Z1M(stpfromMi)+stpZM.*Q(indtmp);
        Z2M(stpfromMi)=Z2M(stpfromMi)+stpZM.^2.*Q(indtmp);
    end
    Q(3,1)=ksfs(1); Q(nxs*nkin^m-2,nxs*nkin^m-8)=ksfs(2); %tweat to be apperiodic
    
    %%%%%detach
    displs=finddispl(ff,ksps,xs(i));
    detfromi=detfrom+(i-1)*nkin^m;
    dettoi=detto+nxs*nkin^m;
    indtmp=sub2ind(size(Q),detfromi,dettoi);
    indtmp=indtmp';
    Q(indtmp(1:3))=kdetfun(displs(1)*ksps(1),kdets(1));   %force dependent
    Q(indtmp(4:6))=kdetfun(displs(2)*ksps(2),kdets(2));
    detZM=xs(i)/2*[-1 -1 -1 1 1 1]';
    Z1M(detfromi)=Z1M(detfromi)+detZM.*Q(indtmp);
    Z2M(detfromi)=Z2M(detfromi)+detZM.^2.*Q(indtmp);
    
    %%%%chemical
    Qblk1=QbABC(displs(1),ksps(1),kabs(1),kbc0s(1),kbas(1),delta1s(1)); %1st motor; displ,ksp,kab,kbc0,kba,delta1
    Qblk2=QbABC(displs(2),ksps(2),kabs(2),kbc0s(2),kbas(2),delta1s(2)); 
    subtmp=((i-1)*nkin^m+1):i*nkin^m;
    Q(subtmp,subtmp)=Q(subtmp,subtmp)+blkdiag(Qblk2,Qblk2,Qblk2)+kron(Qblk1,eye(nkin)); 
end


%%%%one motor attach
inx0=find(xs==0);
Qblk2=QbABC(ff/ksps(2),ksps(2),kabs(2),kbc0s(2),kbas(2),delta1s(2)); %displ,ksp,kab,kbc0,kba,delta1
Qblk1=QbABC(ff/ksps(1),ksps(1),kabs(1),kbc0s(1),kbas(1),delta1s(1));
sub1motor =  (1:6) + nxs*nkin^m;
Q(sub1motor,sub1motor)=Q(sub1motor,sub1motor)+blkdiag(Qblk2,Qblk1);
attto1=inx0+[stptoP stptoM];
ind1motoratt=sub2ind(size(Q),sub1motor,attto1);
Q(ind1motoratt)=katt;
%%%%%detach to zero motor
%Q(nxs*nkin^m+[1 4],end)=kdetfun(ff);

% %%%%zero motor attach, need it???
% Q(end,[1 4]+nxs*nkin^m)=katt;

%%%%%%%%%%%moments
% Z1M=[0;Z1M];
% Z2M=[0;Z2M];
sumrates=sum(Q,2);
M1tau=1./sumrates;
M2tau=2*M1tau.^2;
Z1M=Z1M./sumrates;
Z2M=Z2M./sumrates;
momtauZ=[M1tau M2tau Z1M Z2M]; 

P=Q./repmat(sumrates,1,totst);
P2=P-ones(totst)/(totst);
% Q=Q-diag(dig);
% stdistr=abs(null(Q'));

%%%%%%%%%%find stdistr; two ways
% [stdistr,dummy]=eigs(P',2);
% if dummy(1)<0
%     keyboard
%     stdistr=stdistr(:,2);
% else
%     stdistr=stdistr(:,1);
% end
% if sign(max(stdistr))~=sign(min(stdistr))
%     %keyboard
%     stdistr=stdistr/sum(stdistr);
%     ind=find(stdistr<0);
%     if abs(sum(stdistr(ind)))>1e-10
%         keyboard
%     else
%         stdistr(ind)=0;
%     end
% end
% stdistr=stdistr/sum(stdistr);

IP=eye(totst)-P';
IP=[IP(1:end-1,:); ones(1,totst)];
rhs=[zeros(totst-1,1); 1];
stdistr=IP\rhs;
if sign(max(stdistr))~=sign(min(stdistr))
    keyboard
    stdistr=stdistr/sum(stdistr);
    ind=find(stdistr<0);
    if abs(sum(stdistr(ind)))>1e-10
        keyboard
    else
        stdistr(ind)=0;
    end
end

%%%%%%%%%%%%
meantau=dot(momtauZ(:,1),stdistr);
meanZ=dot(momtauZ(:,3),stdistr);
V(iF)=Lstep*meanZ/meantau;
varZ=findvar(momtauZ(:,3),momtauZ(:,4),meanZ,stdistr,P2,totst);
vartau=findvar(momtauZ(:,1),momtauZ(:,2),meantau,stdistr,P2,totst);
covZtau=findcov(momtauZ(:,1),momtauZ(:,3),meantau,meanZ,stdistr,P2,totst); 
D(iF)=Lstep^2*(meanZ^2*vartau/(meantau)^3+varZ/meantau-2*(meanZ*covZtau/meantau^2))/2;
A=Q-diag(sumrates);
A(:,end)=1;
stdistrcon(:,iF)=A'\[zeros(totst-1,1); 1];
%test=null(Q-diag(sumrates));
%%%%%%%%run length
B=zeros(totst,1);
B(totst-5)=kdetfun(ff,kdets(2));
B(totst-2)=kdetfun(ff,kdets(1));
sumrates2=sumrates+B;
A=Q-diag(sumrates2);
mutau=-(A'\stdistrcon(:,iF))'*ones(totst,1);
ER(iF)=V(iF)*mutau;

end

end