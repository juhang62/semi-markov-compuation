
function out=findvar(M1,M2,meanZ,stdist,P,totxs)
mumu=M1-meanZ;
ImP=eye(totxs)-P;
tmp=ImP\mumu;
covZ=mumu'*diag(stdist)*tmp-mumu'*diag(stdist)*mumu; %last term b/c sum start at 1
out=dot(M2,stdist)-meanZ^2 + 2*covZ;
%out=dot(M2,stdist)-meanZ^2;
end