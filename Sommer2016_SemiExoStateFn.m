function prob=Sommer2016_SemiExoStateFn(n,nprime,K,probstayhome)

% Matlab has a binopdf, but it cannot be used here as it conflicts with arrayfun()
% prob=binopdf(nprime,n+K,probstayhome);

% Following is just a manual verion
% prob= nchoosex * successprob^x * failprob^(n-x)
% prob= nchoosek(n+K,nprime) * probstayhome^nprime * (1-probstayhome)^(n+K-nprime);
% Turns out cannot use nchoosek() with arrayfun() either so,
% Following is a manual verion
% prob= (factorial(n+K)/(factorial(nprime)*factorial(n+K-nprime))) * probstayhome^nprime * (1-probstayhome)^(n+K-nprime);
% Turns out cannot use factorial() with arrayfun() either so,

if nprime>n+K
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


end
