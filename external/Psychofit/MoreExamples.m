% Here are two examples of using the psychofit tools to fit psychometric data. 
% Example 1: data from 0 to 1
% Example 2: data from 0.5 to 1

%% ---------------- Example 1 ----------------
% 
% Fitting data from 0 to 1, using erf_psycho
% Notice that we log the x values before passing them to erf_psycho


%addpath('\\zserver\Code\Psychofit');

% stimulus values
xx = 10.^([-4:0.2:-1]);
nxx = length(xx);

% actual parameters of the observer
threshold = -2;
slope = 1;
gamma = 0.05;

pp = erf_psycho([threshold slope gamma],log10(xx));

% fake experimental data given those parameters
ntrials = 40; 
dd = NaN*zeros(ntrials,nxx);
for istim = 1:nxx
	dd(:,istim) = binornd(1, pp(istim), ntrials, 1);
end

% fit to reconstruct the parameters
pars = mle_fit_psycho([log10(xx); 10*ones(1,nxx); mean(dd)],'erf_psycho')

% graphics

figure; clf
semilogx(xx, mean(dd), 'o');
hold on
semilogx(xx,pp,'--');
semilogx(xx, erf_psycho(pars,log10(xx)), '-');

%% -------------- Example 2 -----------------
%
% Fitting data from .5 to 1, using weibull50
% Notice that we DO NOT log the x values before passing them to erf_psycho
%

addpath('C:/users/matteo/matlab/toolbox5/psychofit');

% stimulus values
xx = 10.^([-4:0.2:-1]);
nxx = length(xx);

% actual parameters of the observer
alpha = 10^(-2.5);
beta = 2;
gamma = 0.05;

pp = weibull50([alpha,beta,gamma],xx);

% fake experimental data given those parameters
ntrials = 40; 
dd = NaN*zeros(ntrials,nxx);
for istim = 1:nxx
	dd(:,istim) = binornd(1, pp(istim), ntrials, 1);
end

% fit to reconstruct the parameters
pars = mle_fit_psycho([xx; 10*ones(1,nxx); mean(dd)],'weibull50')

% graphics

figure; clf
semilogx(xx, mean(dd), 'o');
hold on
semilogx(xx,pp,'--');
semilogx(xx, weibull50(pars,xx), '-');

%% Example 3 -- two different gammas

% make fake data

threshold = -10;
slope = 20;
gamma1 = 0.2;
gamma2 = 0;

xx = -50:10:50;
dd = erf_psycho_2gammas([threshold slope gamma1 gamma2],xx);

% fit it
pars = mle_fit_psycho([xx; repmat(10,[1 11]); dd],'erf_psycho_2gammas',[0 20 0.1 0.1],[-20 10 0 0],[20 40 0.3 0.3],10);

% look at the results
figure; 
plot(xx,dd,'ko','markerfacecolor','k'); hold on
plot(-50:50, erf_psycho_2gammas( pars, -50:50 ));
set(gca,'ylim',[0 1]);


