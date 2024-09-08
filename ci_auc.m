function [ci Serror] = ci_auc(Area, lu, lh)
% confident interval of the AUC
%

alpha = 0.05;   % default
%Area is the AUC (area under ROC)
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(lu-1)*(Q1-Area2)+(lh-1)*(Q2-Area2))/(lu*lh);
Serror=realsqrt(V);
%confidence interval
ci=Area+[-1 1].*(realsqrt(2)*erfcinv(alpha)*Serror);

