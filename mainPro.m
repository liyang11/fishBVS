function [] = mainPro(id)
global X Y W M n a_tau2 invb_tau2 usejeffery

verbose = 0;
updateLambda = 1;
updateAlpha = 0;
usejeffery = 0;
nonspat = 1;
initslambda = 1e5;
B = 1e5; burn = 5e4;
checktie = 0;
getsimudat = 0; nbin = 1e4;%simu = 0;
repl = str2double(num2str(id)); ch0 = repl;
nchain = 10;
ncase = ceil(repl/nchain);
repl = repl - (ncase-1)*nchain;

load('data.mat')

[n,p] = size(X);
pstar = min(p, n-1);
x.kappa = 1-n^(-0.5); % alpha:   0, 1: AIC, log(n)/2: BIC,  or kappa:  0, 1-exp(-1): AIC, 1-n^(-0.5): BIC

if updateAlpha == 1
    kappagap = 1e-2; kappavec = kappagap:kappagap:(1-kappagap);
    EV.kappaconst = log(kappavec) - log(1-(1-kappavec).^(pstar+1));
    EV.kappavec = kappavec; EV.kappagap = kappagap;
end

fprintf('case = %d, chain = %d, sample size = %d, pstar = %d:\n', [ncase, repl, n, pstar])

if nonspat == 1
    M = eye(n); W = zeros(n);
end

% set hyperparamers
meantau2 = .01; vartau2 = 1e2;
a_tau2 = 2+meantau2^2/vartau2;
invb_tau2 = meantau2*(a_tau2-1);
% a_tau2 = 0; invb_tau2 = 0;

meanlambda = .01; varlambda = 1e2; %100
a_lambda = 2+meanlambda^2/varlambda;
invb_lambda = meanlambda*(a_lambda-1);

rng('default'); rng(ch0*12);
% p0s = [2,5,9,12,p-2];
p0s = round(linspace(1,p,nchain)); %[2,5,3,6,pstar-1];
p0 = p0s(repl);
if getsimudat == 1
    p0 = 6; % for simudat
end
x.inds = randsample(1:p, p0);
if all(x.inds~=26)
    x.inds(p0) = 26;
end
x.p = length(x.inds);
x.gamma = 0; x.lambda = initslambda*ones(1, x.p+1);
[x,L0] = getLoglike(x,1);

xall = cell(1,B-burn); ps = zeros(1,B);
births = nan(1,B); deaths = nan(1,B); %swaps = nan(1,B);

% rng(100*repl);
u0 = rand(1,B); u1 = log(rand(1,B)); u2 = rand(1,B);
tic
for b = 1:B
    if verbose == 1
        fprintf('%6d', x.p)
        if ~mod(b,20)
            fprintf('\n')
        end
    end
    
    [x,L0] = getLoglike(x,0);
    
    if (u0(b) <= 0.5 || x.p == 0) && (x.p < pstar) % birth move
        i0 = randsample(setdiff(1:p,x.inds),1);
        x1 = x; x1.inds = [x.inds, i0]; x1.p = length(x1.inds);
        x1.lambda = [x1.lambda, 1/gamrnd(a_lambda, 1/(invb_lambda))];
        [x1,L1] = getLoglike(x1,1);
        logratio = L1 - L0 + log(1-x1.kappa);
        if x.p==0
            logratio = logratio + log(.5);
        end
        MHratio = min(logratio,0); birth = 0;
        if u1(b) <= MHratio
            x = x1; birth = 1; L0 = L1;
        end
        births(b) = birth;
    else % death move
        i0 = randsample(x.p,1);
        x1 = x; x1.inds(i0) = []; x1.p = length(x1.inds); x1.lambda(i0+1) = [];
        [x1,L1] = getLoglike(x1,1);
        logratio = L1 - L0 - log(1-x1.kappa);
        if x.p == pstar
            logratio = logratio + log(2);
        end
        
        MHratio = min(logratio,0); death = 0;
        if u1(b) <= MHratio
            %             if x.inds(optind) == 26
            %                 disp('days_rec is considered out')
            %             end
            x = x1; death = 1; L0 = L1;
        end
        deaths(b) = death;
    end
    
    % update tau2
    Xtild = [ones(n,1),X(:,x.inds)];
    rs = Xtild*x.beta; err = Y - rs;
    if nonspat == 0
        Sigma = M-x.gamma.*W; Lo = chol(Sigma, 'lower');
        L = 0.5*sum((Lo'*err).^2) + 0.5*sum((x.beta).^2./x.lambda');
    else
        L = 0.5*sum((err).^2) + 0.5*sum((x.beta).^2./x.lambda');
    end
    x.tau2 = 1/gamrnd(a_tau2+0.5*(n+x.p+1), 1/(invb_tau2+L));
    
    %update beta
    if nonspat == 1
        Sigma = Xtild'*Xtild + diag(1./x.lambda);
        Mu = Xtild'*Y;
    else
        Sigma = Xtild'*Sigma; Mu = Sigma*Y;
        Sigma = Sigma*Xtild + diag(1./x.lambda);
    end
    Lo = chol(Sigma, 'lower');
    Mu = Lo\Mu;
    Mu = Mu + normrnd(0, sqrt(x.tau2), [length(Mu),1]);
    x.beta = Lo'\Mu;
    
    % update Lambda
    if updateLambda == 1
        for j = 1:(x.p+1)
            x.lambda(j) = 1/gamrnd(a_lambda+1/2, 1/(invb_lambda+sum(x.beta(j).^2)./(2*x.tau2)));
        end
    end
    
    % update Alpha
    if updateAlpha == 1
        loglike = EV.kappaconst + x.p*log(1-EV.kappavec);
        MaxLogLike = max(loglike);
        P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
        U0 = rand(1);
        cump = cumsum([0 P(1:(end-1))]);
        i0 = sum(U0 > cump);
        x.kappa = EV.kappavec(1);
        if i0 > 1
            x.kappa = EV.kappavec(i0-1) + EV.kappagap/P(i0)*(U0-cump(i0));
        end
    end
    
    % update Gamma
    if nonspat == 0
        err = Y - rs;
        tmp = err'*W*err./(2*x.tau2);
        loglike = loglike0 + gammas*tmp;
        MaxLogLike = max(loglike);
        P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
        cump = cumsum([0 P(1:(end-1))]);
        i0 = sum(u2(b) > cump);
        gammanew = gammas(1);
        if(i0>1)
            gammanew = gammas(i0-1) + gaps/P(i0)*(u2(b)-cump(i0));
        end
        x.gamma = gammanew;
    end
    
    if verbose == 1 && b > nbin && ~mod(b,nbin)
        rbirths = births((b-nbin):b); rdeaths = deaths((b-nbin):b); %rswaps = swaps((b-nbin):b);
        fprintf('iter=%d, birth=%.3f, death=%.3f, p=%d, ML = %4.2f, kappa = %4.2f, tau2 = %4.2f, lambda = %4.2f\n',...
            [b, mean(rbirths(~isnan(rbirths))), mean(rdeaths(~isnan(rdeaths))), x.p, L0, x.kappa, x.tau2, mean(x.lambda)]);
        %swap=%.3f,   mean(rswaps(~isnan(rswaps))),
    end
    
    ps(b) = x.p;
    if checktie == 1 && b>21
        if(all(ps((b-20):b) > 10))
            disp('here')
        end
    end
    
    % record results
    if(b > burn)
        xall{b-burn} = x;
        %Zmat(b-burn,:) = Y(zind)';
    end
end

CPUtime = toc; CPUtime = CPUtime/60;
fprintf('%d iterations are done with elapsed time %.2f minutes.\n', B, CPUtime)
nam = strcat('MSCopt',num2str(ncase), '_', num2str(repl),'.mat');
save(nam,'xall','CPUtime','births','deaths','ps','nonspat')
end

function [x,L] = getLoglike(x, propose)
% for flat prior
global X Y n W M a_tau2 invb_tau2 usejeffery
Xtild = [ones(n,1),X(:,x.inds)];
Sigma = M-x.gamma.*W; L = Y'*Sigma*Y;
Sigma = Xtild'*Sigma; Mu = Sigma*Y;
Sigma = Sigma*Xtild + diag(1./x.lambda);
Lo = chol(Sigma, 'lower');
Mu = Lo\Mu; L = (L - sum(Mu.^2))/2;
n1 = 0.5*n;
if propose == 1
    x.tau2 = 1/gamrnd(a_tau2+n1, 1/(invb_tau2+L));
    Mu = Mu + normrnd(0, sqrt(x.tau2), [length(Mu),1]);
    x.beta = Lo'\Mu;
end
if usejeffery == 0
    L = - 0.5*sum(log(x.lambda)) ...
        - sum(log(diag(Lo))) - (a_tau2+n1)*log(1+L/invb_tau2);
elseif usejeffery == 1
    L = - 0.5*sum(log(x.lambda)) ...
        - sum(log(diag(Lo))) - (a_tau2+n1)*log(L);
end
end
