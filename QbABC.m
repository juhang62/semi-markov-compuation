function A=QbABC(displ,ksp,kab,kbc0,kba,delta1) %transition that leads no change in config; kinABC; ff +hindering
global kbt
%%%%need cutoff!!!!!!!!!!!
%fdep = @(f) min(exp(-f*delta1/kbt),1e5);%*double(displ+1<=3); make not too big
fdep = @(f) exp(-f*delta1/kbt);%*double(displ+1<=3); 
kbc=kbc0*fdep(ksp*displ);%displ*ksp);%displ*ksp);
A=[0 kab 0; kba 0 kbc;0 0 0];

end