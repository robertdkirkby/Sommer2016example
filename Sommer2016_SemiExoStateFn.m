function prob=Sommer2016_SemiExoStateFn(n,nprime,K,probstayhome,maxnumchildren)
% Transition probability for child count nprime conditional on n and K.
% Inputs are scalars; output is the scalar transition probability.

% Matlab has a binopdf, but it cannot be used here as it conflicts with arrayfun()
% prob=binopdf(nprime,n+K,probstayhome);

% Following is just a manual verion
% prob= nchoosex * successprob^x * failprob^(n-x)
% prob= nchoosek(n+K,nprime) * probstayhome^nprime * (1-probstayhome)^(n+K-nprime);
% Turns out cannot use nchoosek() with arrayfun() either so,
% Following is a manual verion
% prob= (factorial(n+K)/(factorial(nprime)*factorial(n+K-nprime))) * probstayhome^nprime * (1-probstayhome)^(n+K-nprime);
% Turns out cannot use factorial() with arrayfun() either so,
%
% VFI Toolkit evaluates this transition for every value on the decision
% grid, including K=1 when n is already at maxnumchildren. The return
% function rules out that choice, but the transition matrix still needs to
% be a valid stochastic matrix during setup.
if n+K>maxnumchildren
    if nprime==maxnumchildren
        prob=1;
    else
        prob=0;
    end
elseif nprime>n+K
    prob=0;
else
    % factorial(nprime)
    f1=1;
    for ii=1:nprime
        f1=f1*ii;
    end
    % factorial(n+K)
    f2=1;
    for ii=1:(n+K)
        f2=f2*ii;
    end
    % factorial(n+K-nprime)
    f3=1;
    for ii=1:(n+K-nprime)
        f3=f3*ii;
    end
    % (f2/(f1*f3)) gives nchoosek(n+K,nprime)
    prob= (f2/(f1*f3)) * probstayhome^nprime * (1-probstayhome)^(n+K-nprime);

    % Finally, something arrayfun() is happy with :)
end


end %end function
