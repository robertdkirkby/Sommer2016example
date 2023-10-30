function F=Sommer2016_ReturnFn(l,x,K,aprime,a,n,epsilon,f,upsilon,agej,Jr,h_j,r,mew,psi1,psi2,theta,gamma,zeta,kappa,qbar,maxnumchildren,pension,earningsshifter)

F=-Inf;

% Note: the constraints (equation 2 and 3 of Sommer (2016)), are not dependent on f
if agej<Jr % working age
    c=a+earningsshifter*exp(h_j+epsilon+upsilon)*(1-l)-x-aprime/(1+r); 
    % Note: exp(h_j+epsilon+upsilon) is what Sommer (2016) denotes w
    % Comment: this aprime/(1+r), instead of just (1+r)*a comes because of how Sommer (2016) solves the 
    % problem but makes no sense from the perspective of how VFI Toolkit nor makes much sense in terms 
    % of the economics. [Sommer (2016) uses the common 'wealth-on-hand' reformulation.]
else % agej>=Jr, so retired
    c=a+pension-aprime/(1+r);
end

% production fn for childrens quality
if n>0
    q=(mew*(x/(n^psi1))^theta + (1-mew)*(l/(n^psi2))^theta)^(1/theta);
else
    q=0; % Not used for anything anyway
end
% q=1; %  DEBUG

if c>0
    % u(c,n,q)
    U_c=(c^(1-gamma))/(1-gamma);
    if n>0
        U_nq=zeta*((n*q)^(1-kappa))/(1-kappa);
    else
        U_nq=0;
    end
    F=U_c+U_nq;
end

% % DEBUG
% if f==1 && K==0 && n<maxnumchildren
%     F=-Inf;
% end

% % DEBUG
% if aprime>0.1
%     F=F+10;
% end

if n>0 && q<qbar
    F=-Inf; % must have q>qbar if there are children (qbar is the "lower bound on children's consumption")
end

if f==0 && K==1 % infertile, so cannot have children
    F=-Inf; % Cannot have children if infertile
end

if K==1 && n==maxnumchildren
    F=-Inf; % You are not allowed to try to leave the top of the grid on number of children 
    % This is just needed because of how the transition probabilities on n
    % are determined (they would otherwise fail to sum to one).
    % Not clear from paper how Sommer (2016) deals with this as I imagine she must have had a maximum n in the codes?
end

end