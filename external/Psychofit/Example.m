%% simulate psychometric data and fit them with Weibull function

% stimulus values
nx = 8;
xx = 10.^(linspace(-4,-1,nx));

%% actual parameters of the observer

alpha = 10^(-2.5);
beta = 2;
gamma = 0.2;

pp = weibull50([alpha,beta,gamma],xx);

%%  fake experimental data given those parameters

ntrials = 200; 
dd = NaN*zeros(ntrials,nx);
for istim = 1:nx
	dd(:,istim) = binornd(1, pp(istim), ntrials, 1);
end

%% look at the data

% estimate error bars based on binomial distribution
p_est = zeros(nx,1);
p_conf = zeros(nx,2);
for ix = 1:nx
    [p_est(ix), p_conf(ix,:)] = binofit(nnz(dd(:,ix)),ntrials,0.05);
end

figure; clf
for ix = 1:nx
    semilogx( xx(ix)*[1 1], p_conf(ix,:),'k-'); hold on
end
semilogx(xx, mean(dd), 'o', 'markerfacecolor','k');

manyxx = 10.^linspace(-4,-1);
semilogx(manyxx, weibull50([alpha,beta,gamma],manyxx),'k--');

%% fit to reconstruct the parameters

FractionCorrect = mean(dd);
FractionCorrect(7) = NaN; % this condition was "not run"
pars = mle_fit_psycho([xx; 10*ones(1,nx); FractionCorrect],'weibull50');

%% graphics

figure; clf
for ix = 1:nx
    semilogx( xx(ix)*[1 1], p_conf(ix,:),'k-'); hold on
end
semilogx(xx, mean(dd), 'o', 'markerfacecolor','k'); 

manyxx = 10.^linspace(-4,-1);
semilogx(manyxx, weibull50([alpha,beta,gamma],manyxx),'k--');
semilogx(manyxx, weibull50(pars,manyxx), 'k-', 'linewidth', 2);

%% confidence intervals for parameters and curves

nboots = 200;
bootpars = zeros(nboots,3);
for iboot = 1:nboots
    ii = ceil(ntrials*rand(ntrials,1));
    bootpars(iboot,:) = ...
        mle_fit_psycho([xx; 10*ones(1,nx); mean(dd(ii,:))],'weibull50');
end

% this would give you a confidence interval
alphas = bootpars(:,1);
figure; hist(log10(alphas))
title(sprintf('alpha = %2.2f +- %2.2f', mean(log10(alphas)), std(log10(alphas))));

DensityMatrix = zeros( 40, length(xx) );
for iboot = 1:nboots
    mm =  ceil(40*weibull50(bootpars(iboot,:),xx));
    for ix = 1:length(xx)
        DensityMatrix( mm(ix), ix) = DensityMatrix( mm(ix), ix) + 1;
    end
end
DensityMatrix = DensityMatrix/nboots;
figure; clf;
imagesc(log10(xx),linspace(0,1,40),DensityMatrix);
set(gca,'ydir','normal');
colormap bone
hold on
plot( log10(xx),pp, 'r' );


