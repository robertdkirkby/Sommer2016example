% Replication of Sommer (2016) - Fertility Choice in a life cycle model with idiosyncratic uninsurable earnings risk
% WARNING: This will not actually give the results of Sommer (2016) as the
% parameter h_j (deterministic component of earnings that depends on age)
% has been lost to the sands of time.
% WARNING: The grids are small, especially the one on assets. You will need
% to increase the grid size if you want accurate results. This code is
% intended more as illustration.

% NOTE: Around line 110 I fit an exponential function, you will need the Matlab 'Curve Fitting Toolbox' to run this.

% Note: I do not have the original values of h_j, the (age-dependent) deterministic
% component of ('hourly') wages. I use a similar calibrated term from
% Kaplan (2012) but of course this is not the same thing and so is not
% going to deliver the exact same answers. (Sommer, 2016, calls this h_t.
% She explains how to calculate it from the data but does not report the
% acutal values.)

% Note: Sommer (2016) reports 'psi' as the 'preference scale' relating to
% children, but in the utility function (bottom page 33) it is denote zeta.
% I follow the utility function and call it zeta here.

% Comment: because infertile is an absorbing state in principle what you
% should do is first solve the infertile problem, and then solve the fertile
% problem. VFI Toolkit cannot handle this level of subtlety and so it is
% just going to solve the fertile and infertile and the same time and put
% lots of zeros for the probability of infertile-->fertile when calculating
% expections. This is a little wasteful computationally but we get away with it :)

% Comment: Sommer (2016) Appendix seems to say she does pure discretization (which
% is what VFI Toolkit does). She says she uses 7 grid points for epsilon
% and 2 for upsilon; I switch to 3 for upsilon as that is minimum for Farmer-Toda discretization. 
% She simulated the agents distribution based on 10,000 household simulations. 
% She does not mention number of grid points used for assets, spending-on-child, and leisure.

n_d=[11,11,2]; % leisure, spending-on-child, child (l,x,K in original notation)
n_a=101; % assets (a in original notation)
n_z=[7,2]; % persistent component of wage, fertile (epsilon, f in original notation)
n_e=3; % transitory component of wage (upsilon in original notation)
n_semiz=6; % number of children (n in original notation) [Note: max number of children will be n_semiz-1 as first point is zero]

N_j=63; % From age 18 to 80
Params.agejshifter=17; % Only used when creating graphs

Params.earningsshifter=1; % Just used to see how changing h_t up/down changes the results (I introduced this just because the actual values of h_t are not available)

%% Parameters
% Discount factor
Params.r=0.07; % 0.04
Params.beta=1/(1+Params.r); % Discount factor

% Preferences
Params.gamma=1.5; % risk aversion

% Child preferences
Params.kappa=0.14; % preference curvature
Params.zeta=3.5; % preference scale
Params.mew=0.35; % production share
Params.theta=0.7; % 1/(1-theta) is the elasticity of substitution in production
Params.qbar=0.34; % lower bound on children's consumption
Params.psi1=0.91; % Household economies to money input to production
Params.psi2=0.54; % Household economies to time input to production
Params.probstayhome=0.98; % Sommer (2016) calls this p.
                         % Table 2 of Sommer (2016) says that 'probability that a child stays home is (1-p)'. But this is
                         % in conflict with her equation 4, which has p being the probability the child stays home, and
                         % 1-p being the probability the child leaves home.  Looking at how equation 4 would evolve if the 
                         % probability of staying home was 0.02 (almost all households instanly drop to 0 children) I am 
                         % confident this was a typo in the table 4 description on what the coefficient represents.
Params.maxnumchildren=n_semiz-1; % Used to avoid people trying to leave the top of the grid on number of children

% Retirement
Params.Jr=66-Params.agejshifter; % Retirement period, Sommer (2016) calls this R; on pg 30 Sommer says retirement is age >65
Params.pension=1; % Just an initial guess (I replace this below)
Params.pensionreplacementrate=0.4; % This is the target for pension

%
Params.agej=1:1:N_j;

%% Wage process
% w follows an exogenous process with three components: deterministic, persistent, and transitory
%
% The deterministic is a age-dependent polynomial (Sommer (2016) explains
% how she gets this, but since the numbers are not reported I was lazy and
% just used the one from Kaplan (2012) which serves a similar purpose.)
% The persistent and transitory are easy enough.
% Because the mean of the Kaplan (2012) values do not fit, I renormalize so
% that mean earnings are equal to 2/3rds of 2004 GDP per capita (2/3 is
% roughly the labor share). This is all rather lazy of me but I don't feel
% like going and doing all the CPS data (with associated complicated
% weights) just to get a deterministic age profile for an example code.
%
% h_j, the age-dependent deterministic component of wages
% Params.h_j= % Sommer (2016) calls this h_t, but I rename to h_j to fit VFI Toolkit convention of using j for age/period
%   Sommer (2016), pg 33: "The average age-profile for wages, h(t), is calculated from the 2004 CPS by 
%   dividing the family labor income, defined as a sum of yearly earnings of both spouses in husband–wife 
%   families, by the sum of total hours worked by the couple. The average age of the couple is taken to 
%   represent the age of the household. The profile is smoothed using a cubic polynomial in age."
% This would be a pain as you need to use the CPS weights, so instead I just take the values used by 
% Kaplan (2012). Since Kaplan (2012) numbers start at age 29 I just use a
% linear trend on 29-35 to backfill to age 18. This is not a good idea, but
% I just want to get this done. I also double the age 64 number as Kaplan
% has retire at 65, while Sommer has retire >65.

% Kaplan (2012) deterministic function of age: ages 29 to 64 inclusive
Params.h_j=[2.576039, 2.620736, 2.663742, 2.705041, 2.744625, 2.782487, 2.818624, 2.853036, 2.885727, 2.916703, 2.945975, 2.973556, 2.999462, 3.023714, 3.046335, 3.067351, 3.086793, 3.104694, 3.12109, 3.136021, 3.149529, 3.161663, 3.172471, 3.182006, 3.190324, 3.197486, 3.203554, 3.208595, 3.212677, 3.215874, 3.218262, 3.219921, 3.220932, 3.221383, 3.221363, 3.220963]; % Contents of 'kappasmooth.txt' from Kaplan (2012) files
% Fill in first part
% slopeinage=(Params.h_j(5)-Params.h_j(1))/4;
slopeinage=0;
Params.h_j=[slopeinage*(-11:1:-1)+Params.h_j(1),Params.h_j,Params.h_j(end)];
% Note, these values are just the working age ones, need to fill in retirement as zeros
Params.h_j=[Params.h_j,zeros(1,N_j-Params.Jr+1)];
% I use something to shift the exp(h_j) to be tiny as this seems to fit what Sommer (2016) used
% Params.earningsshifter=0.00002*3/20; % The 0.00002 is roughly what we need c to be, the /20 is rough size of what exp(h_j) would otherwise be

% Update pension guess (remains a guess)
Params.pension=Params.pensionreplacementrate*0.5*exp(Params.h_j(1)); % the 0.5 was to allow for some leisure



% The persistent and transitory components are just following Sommer
% (2016), except I use Farmer-Today method to discretize (I expect she used
% Tauchen, Farmer-Toda is better but anyway did not exist in 2016).
% epsilon, the persistent AR(1) component
Params.rho_w=0.56;
Params.sigma_epsilon=0.21;
[epsilon_grid,pi_epsilon] = discretizeAR1_FarmerToda(0,Params.rho_w,Params.sigma_epsilon,n_z(1));
% upsilon, the transitory iid component
Params.sigma_upsilon=0.17;
[upsilon_grid,pi_upsilon] = discretizeAR1_FarmerToda(0,0,Params.sigma_upsilon,n_e);
pi_upsilon=pi_upsilon(1,:)'; % iid, so just set as a column vector

%% Grids
l_grid=linspace(0,1,n_d(1))'; % labor supply as fraction of time
x_grid=linspace(0,1,n_d(2))'; % money spent on child quality
K_grid=[0;1]; % 0=not have a baby, 1=have a baby

a_grid=200*linspace(0,1,n_a).^3'; % the power of 3 means points are denser at low asset levels which is where value function has most curvature and so tends to be good for accuracy

% In principle, the grid on fertility is just
% f_grid=[0,1]; % 0=infertile, 1=fertile
% But because the transition probabilities for f are age-dependent, VFI Toolkit
% requires the grid to also be set up as age-dependent (even though in this
% instance it is not actually changing with age)
f_grid_J=[0;1].*ones(1,N_j); % 0=infertile, 1=fertile
% Sommer (2016) does not report the actual values of the transition
% probabilities. She does clearly explain how they were calculated so that
% is now repeated here. First, we get the data points from Trussell &
% Wilson (1985, https://doi.org/10.1080/0032472031000141486 ) that she
% shows in Figure 1 on pg 32. These are
TW1985values=[7.0,13.1,23.1,34.5,57.6,95.2]';
TW1985ages=[20,25,30,35,40,45]';
% Reading Trussell & Wilson (1985) with an eye to finding something that
% looks like Fig 1 of Sommer (2016) it looks to me like she uses Figure 2,
% the line corresponding to women who married at age 15-19. This curve
% appears to be identical to that in Fig 1 of Sommer (2016). Footnote of
% the TW1985 figure says that the source for this is Tables 4 and 8 of
% their paper. In Table 8 of TW1985, this seems to correspond to the first
% row, so I here use those numbers. Sommer (2016) has first point in her
% figure as age 20, this appears to correspond to the 20-24 age grouping of
% TW1985 (so skip past first entry of row 1 with corresponds to age 15-19).
% Next step is Sommer (2016) reports fitting an exponentional function of
% age to this.
expf = fit(TW1985ages,TW1985values,'exp1'); % NOTE: You need the Matlab 'Curve Fitting Toolbox' to run this
% Let's do a version of Figure 1 of Sommer (2016) to make sure things look okay so far
fig1=figure(1);
plot(expf,TW1985ages,TW1985values)
ylabel('Fraction of Couples Permanently Infertile')
xlabel('Age of Woman')
title('Figure 1 of Sommer (2016)')
% Sommer (2016) then sets the transition probabilities for f, which she
% denotes p^I_t to match this: "The probabilities, p^I_t, are derived so
% that the fraction of permanently infertile households of any given age in
% the model matches exactly the corresponding fraction in the data. In the
% data about 97 percent of all couples are infertile at age 45. In the
% model, the cumulative probability that a household is permanently
% infertile at age 45 is set to 1."
% So first, use our fitted exponential to get the fraction infertile for each age.
frac_infertile=expf(Params.agejshifter+1:45)/100; % Note: now need actual units to be fraction
f_initialdist=[frac_infertile(1); 1-frac_infertile(1)]; % Initial distribution in first period
pi_f_J=zeros(2,2,N_j);
% Once we reach age 45 everyone is infertile
for jj=(45-Params.agejshifter):N_j
    pi_f_J(:,:,jj)=[1,0;1,0]; % probability 1 of shifting to infertile
end
% Before that, infertile is absorbing and need probability
% fertile-infertile to deliver the next age frac_infertile
for jj=1:(45-Params.agejshifter-1)
    % first row: infertile is absorbing state
    % second row: next period infertile = this period infertile + pi_fertiletoinfertile*this period fertile
    pi_fertiletoinfertile=(frac_infertile(jj+1)-frac_infertile(jj))/(1-frac_infertile(jj));
    pi_f_J(:,:,jj)=[1,0; pi_fertiletoinfertile, 1-pi_fertiletoinfertile]; 
end
% To double-check, simulate the fertilite process and make sure we get the correct fraction infertile
pi_f_dist=zeros(2,N_j);
pi_f_dist(:,1)=f_initialdist;
for jj=1:(N_j-1)
    pi_f_dist(:,jj+1)=pi_f_J(:,:,jj)'*pi_f_dist(:,jj);
end
checka=pi_f_dist(1,1:45-Params.agejshifter)-frac_infertile'; % This should all be zeros, it is
checkb=pi_f_dist(1,45-Params.agejshifter+1:end); % This should all be ones, it is
% Sucess :)

% Now put the grids together in right form for VFI Toolkit
d_grid=[l_grid; x_grid; K_grid]; % leisure, spending-on-child, child (l,x,K in original notation)
% a_grid
z_grid_J=[epsilon_grid.*ones(1,N_j); f_grid_J]; % Because f transitions depend on j, have to make all grids and transition matrices depend on j
pi_z_J=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    % Recall: z=(epsilon,f)
    pi_z_J(:,:,jj)=kron(pi_f_J(:,:,jj),pi_epsilon);
    % Use kron() in reverse order on the shocks that make up z
end
% Need to set up an age-dependent grid on z by putting it into vfoptions
% and simoptions, and then putting a placeholder into its place
z_grid=z_grid_J(:,1); % is just a placeholder
vfoptions.z_grid_J=z_grid_J;
pi_z=pi_z_J(:,:,1); % is just a placeholder
vfoptions.pi_z_J=pi_z_J;
% e variables have to go into vfoptions
vfoptions.n_e=n_e;
vfoptions.e_grid=upsilon_grid;
vfoptions.pi_e=pi_upsilon;
% Also have to put these into simoptions
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;
simoptions.z_grid_J=vfoptions.z_grid_J;
simoptions.pi_z_J=vfoptions.pi_z_J;

%% Children, semi-exogenous state
% Just need to set up the transitions on n, which evolve as Binominal(n+K,p) [equation 4 of Sommer (2016)]

n_grid=(0:1:(n_semiz-1))'; % first grid point is zero children

% Set up the semi-exogenous state
vfoptions.n_semiz=n_semiz;
vfoptions.semiz_grid=n_grid;
% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(n,nprime,K,probstayhome) Sommer2016_SemiExoStateFn(n,nprime,K,probstayhome);
% It is hardcoded that only the 'last' decision variable can influence the transition probabilities of the semi-exogenous states
% The semi-exogenous states must be included in the return fn, fns to evaluate, etc., as the last 'z' variables


%% Return function
DiscountFactorParamNames={'beta'};

% order must be: (d,aprime,a,z,semiz,e,...)
ReturnFn=@(l,x,K,aprime,a,epsilon,f,n,upsilon,agej,Jr,h_j,r,mew,psi1,psi2,theta,gamma,zeta,kappa,qbar,maxnumchildren,pension,earningsshifter) ...
    Sommer2016_ReturnFn(l,x,K,aprime,a,epsilon,f,n,upsilon,agej,Jr,h_j,r,mew,psi1,psi2,theta,gamma,zeta,kappa,qbar,maxnumchildren,pension,earningsshifter);

%% Solve for the value function and policy fn
vfoptions.verbose=1;
vfoptions.lowmemory=1
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z, N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
timevf=toc


%% Initial distribution
% f_initialdist was defined above
% start with n=0; (Sommer pg 30, n_18=0)
% start with epsilon=0; (Sommer pg 30, epsilon_1=0)
% start with upsilon=0; % (Sommer does not appear to specify, but I am
% Sommer (2016) Appendix says she used 2 points for epsilon, so I assume that just start on the stationary dist for epsilon.
jequaloneDist=zeros([n_a,n_z,n_semiz,n_e],'gpuArray');
jequaloneDist(1,ceil(n_z(1)/2),:,1,:)=shiftdim(f_initialdist.*shiftdim(vfoptions.pi_e,-2),-3); % This would be easier to set up as a loop, but I can just do it as a high-dimension product

%% Population weights
% There is no death (until final period) and no population growth, so just use equal weights.
Params.mewj=ones(1,N_j)/sum(N_j);
AgeWeightsParamNames={'mewj'};

%% Solve for the stationary dist
% We also need to tell simoptions about the semi-exogenous states
simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
% Need to also tell simoptions about the semi-exogenous shocks
% Because evaluating pi_semiz_J requires the d_grid we also have to provide
simoptions.d_grid=d_grid;

simoptions.verbose=1;
simoptions.parallel=5;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% FnsToEvaluate are how we say what we want to graph the life-cycles of
% Like with return function, we have to include (h,aprime,a,z) as first
% inputs, then just any relevant parameters.
FnsToEvaluate.leisure=@(l,x,K,aprime,a,epsilon,f,n,upsilon) l; % l is fraction of time for leisure
FnsToEvaluate.moneyonkids=@(l,x,K,aprime,a,epsilon,f,n,upsilon) x; % x is money spend on children
FnsToEvaluate.havechild=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) K; % K is decision to have a child
FnsToEvaluate.assets=@(l,x,K,aprime,a,epsilon,f,n,upsilon) a; % a is the current asset holdings
FnsToEvaluate.fractiontimeworked=@(l,x,K,aprime,a,epsilon,f,n,upsilon) 1-l; % l is fraction of time for leisure
FnsToEvaluate.earnings=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) exp(h_j+epsilon+upsilon)*(1-l); % w*kappa_j*z*(1-l) is the labor earnings
FnsToEvaluate.nchildren=@(l,x,K,aprime,a,epsilon,f,n,upsilon) n; % n is number of children
FnsToEvaluate.infertility=@(l,x,K,aprime,a,epsilon,f,n,upsilon) (1-f); % state that determines whether it is possible to have children
FnsToEvaluate.childquality=@(l,x,K,aprime,a,epsilon,f,n,upsilon,mew,psi1,theta,psi2) (n>0)*((mew*(x/(n^psi1))^theta + (1-mew)*(l/(n^psi2))^theta)^(1/theta)); % q, child quality
% notice that we have called these fractiontimeworked, earnings, assets and nchildren

FnsToEvaluate.condlearnings=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) exp(h_j+epsilon+upsilon); % w*kappa_j*z*(1-l) is the labor earnings
FnsToEvaluate.h_j=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) h_j; % w*kappa_j*z*(1-l) is the labor earnings
FnsToEvaluate.epsilon=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) epsilon; % w*kappa_j*z*(1-l) is the labor earnings
FnsToEvaluate.upsilon=@(l,x,K,aprime,a,epsilon,f,n,upsilon,h_j) upsilon; % w*kappa_j*z*(1-l) is the labor earnings


%% Calculate the life-cycle profiles
simoptions.parallel=2;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

% We can do a version of Figure 2 of Sommmer (2016)
fig2=figure(2);
plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.nchildren.Mean)
xlim([18,45]) % just so it matches Sommer (2016)
ylabel('Age')
xlabel('Average Number of Children Born') % really it is just number of children in the household (not sure if Sommer (2016) actually plots the number born or just the number; number born would be easy enough, just set up a FnsToEvaluate that returns K, get the life-cycle profile for the mean of K, then to the cumulative sum of age-conditional-mean-K over age, then plot)
title('Figure 2 of Sommer (2016)')

fig3=figure(3);
subplot(4,2,1); plot(AgeConditionalStats.fractiontimeworked.Mean)
title('Fraction time worked (1-l)')
subplot(4,2,2); plot(AgeConditionalStats.moneyonkids.Mean)
title('Money on kids (x)')
subplot(4,2,3); plot(AgeConditionalStats.havechild.Mean)
title('Have a child (K)')
subplot(4,2,4); plot(AgeConditionalStats.assets.Mean)
title('Assets (a)')
subplot(4,2,5); plot(AgeConditionalStats.earnings.Mean)
title('Earnings')
subplot(4,2,6); plot(AgeConditionalStats.infertility.Mean)
title('Infertility (1-f)') % Note, this is just plotting Figure 1 of Sommer (2016) again
subplot(4,2,7); plot(AgeConditionalStats.childquality.Mean)
title('Child quality (q)')
subplot(4,2,8); plot(AgeConditionalStats.nchildren.Mean)
title('Number of children (n)')

fig4=figure(4);
subplot(3,2,1); plot(AgeConditionalStats.earnings.Mean)
title('Earnings')
subplot(3,2,2); plot(AgeConditionalStats.fractiontimeworked.Mean)
title('Fraction time worked (1-l)')
subplot(3,2,3); plot(AgeConditionalStats.earnings.Mean)
title('condlearnings')
subplot(3,2,4); plot(AgeConditionalStats.h_j.Mean)
title('h_j')
subplot(3,2,4); plot(AgeConditionalStats.epsilon.Mean)
title('epsilon')
subplot(3,2,4); plot(AgeConditionalStats.upsilon.Mean)
title('upsilon')


%% Set up the pension target as a general equilbrium condition and then solve to get that correct
% Since the calibration is anyway a mess I figured there was no point doing this



% %% Utility of one child versus none
% % % n=1;
% % % q=Params.qbar; % the min
% % % U_nq=Params.zeta*((n*q)^(1-Params.kappa))/(1-Params.kappa);
% % % % note that U_nq with 0 children is just 0
% % % % So one child gives roughly U_nq=459
% % % 
% % % % Imagine we pay cash for kid, so
% % % %     q=(mew*(x/(n^psi1))^theta + (1-mew)*(l/(n^psi2))^theta)^(1/theta);
% % % % with l=0 gives
% % % % q=mew*(x/(n^psi1))
% % % % and n=1 so just q=mew*x
% % % % so x=0.97 (for min q=0.34)
% % % % During working age, x is one for one with consumption
% % % % so we need to find c such that change in c of 0.97 gives a change in U(c)
% % % % of 459
% % % % so (c+0.97)^(1-1.5)-c^(1-1.5)=-0.5*459
% % % testfn=@(c) (c+0.97).^(1-1.5)-c.^(1-1.5)
% % % testfn2=@(c) abs((c+0.97).^(1-1.5)-c.^(1-1.5) +0.5*459)
% % % 
% % % % The magnitude of c such that testfn2 evaluates to zero is around 0.00002
% % % testfn2(0.00001:0.00001:0.00007)
% % % 
% % % 
% % % 
% 
% l_val=l_grid(1);
% x_val=x_grid(1);
% 
% aprime_val=a_grid(1);
% a_val=a_grid(1);
% 
% n_val=n_grid(1);
% 
% epsilon_val=epsilon_grid(4);
% upsilon_val=upsilon_grid(2);
% f_val=1; % fertile
% 
% ageind=63;
% 
% K_val=K_grid(1);
% Sommer2016_ReturnFn(l_val,x_val,K_val,aprime_val,a_val,epsilon_val,f_val,n_val,upsilon_val,Params.agej(ageind),Params.Jr,Params.h_j,Params.r,Params.mew,Params.psi1,Params.psi2,Params.theta,Params.gamma,Params.zeta,Params.kappa,Params.qbar,Params.maxnumchildren,Params.pension,Params.earningsshifter)
% K_val=K_grid(2);
% Sommer2016_ReturnFn(l_val,x_val,K_val,aprime_val,a_val,epsilon_val,f_val,n_val,upsilon_val,Params.agej,Params.Jr,Params.h_j,Params.r,Params.mew,Params.psi1,Params.psi2,Params.theta,Params.gamma,Params.zeta,Params.kappa,Params.qbar,Params.maxnumchildren,Params.pension,Params.earningsshifter)
% 
% temp2=Policy(3,:,4,1,1,2,1);
% [min(temp2(:)),max(temp2(:))] % have child
% temp2=Policy(3,:,4,1,6,2,1);
% [min(temp2(:)),max(temp2(:))]
% 
% temp2=Policy(2,:,4,1,1,2,1);
% [min(temp2(:)),max(temp2(:))] % spending on child
% temp2=Policy(2,:,4,1,6,2,1);
% [min(temp2(:)),max(temp2(:))]
% 
% temp2=Policy(4,:,4,1,1,2,49);
% [min(temp2(:)),max(temp2(:))] % assets
% temp2=Policy(4,:,4,1,6,2,49);
% [min(temp2(:)),max(temp2(:))]
% 
% 
% % no kids
% subplot(4,1,1); plot(1:1:n_a,V(:,4,1,1,2,63),1:1:n_a,V(:,4,2,1,2,63))
% % one kid
% subplot(4,1,2); plot(1:1:n_a,V(:,4,1,2,2,63),1:1:n_a,V(:,4,2,2,2,63))
% % two kids
% subplot(4,1,3); plot(1:1:n_a,V(:,4,1,3,2,63),1:1:n_a,V(:,4,2,3,2,63))
% % three kids
% subplot(4,1,4); plot(1:1:n_a,V(:,4,1,4,2,63),1:1:n_a,V(:,4,2,4,2,63))
% 
% 
% V(1,4,2,1,2,63)


