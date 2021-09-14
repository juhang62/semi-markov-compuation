
function out=findcov(M1v1,M1v2,meanv1,meanv2,stdist,P,realtotxs)
mumuvar1=M1v1-meanv1;
mumuvar2=M1v2-meanv2;
ImP=eye(realtotxs)-P;
covii=dot(M1v1.*M1v2,stdist)-meanv1*meanv2; %Z and tau indepdent given S
cov0=mumuvar1'*diag(stdist)*mumuvar2;
cov1=mumuvar1'*diag(stdist)*(ImP\mumuvar2);
cov2=mumuvar2'*diag(stdist)*(ImP\mumuvar1);
out=covii-2*cov0+cov1+cov2;
end