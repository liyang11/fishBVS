function [] = sumPro()
simu = 0; collectit = 1;
load('data.mat')
Y0 = Y;

for ncase = 1:2
    fprintf('case = %d\n', ncase)
    
    if collectit == 1
        [n,p] = size(X);
        if simu == 1
            load('simudat.mat')
        end
        
        chs = [1:5];
        nch = length(chs);
        niter = 5e4; burn = 0; thin = 5; realburn = 5e4; %5e4
        nsample = (niter-burn)/thin;    tot = nch*nsample;
        allps = zeros(nsample,nch); allgammas = zeros(nsample,nch); alltau2s = zeros(nsample,nch);
        Allps = zeros(realburn+burn+niter, nch);
        allbirths = zeros(nsample,nch); alldeaths = zeros(nsample,nch);
        indmat = zeros(tot,p); betamat = zeros(tot,p); beta0mat = zeros(tot,1);
        for ch = 1:nch
            load(strcat('MSCopt',num2str(ncase),'_',num2str(chs(ch)),'.mat'))
            Allps(:,ch) = ps;
            allps(:,ch) = ps(realburn+burn+(1:nsample).*thin);
            allbirths(:,ch) = births(realburn+burn+(1:nsample).*thin);
            alldeaths(:,ch) = deaths(realburn+burn+(1:nsample).*thin);
            for i = 1:nsample
                indmat((ch-1)*nsample+i,xall{burn+i*thin}.inds) = 1;
                beta0mat((ch-1)*nsample+i) = xall{burn+i*thin}.beta(1);
                betamat((ch-1)*nsample+i,xall{burn+i*thin}.inds) = xall{burn+i*thin}.beta(2:end)';
                allgammas(i,ch) = xall{burn+i*thin}.gamma;
                alltau2s(i,ch) = xall{burn+i*thin}.tau2;
            end
        end
        
        save(strcat('sumout',num2str(ncase),'.mat'),'allps','tot','allbirths','alldeaths','indmat','betamat','allgammas','alltau2s')
    end
    
    load(strcat('sumout',num2str(ncase),'.mat'))
    
    % plot(Allps)
    plot(allps) %(:,[1,2,4,5,6])
    % plot(allgammas)
    ps = reshape(allps, [1,tot]);
    tabulate(ps)
    
    mbirth = nanmean(reshape(allbirths, [1,tot]));
    mdeath = nanmean(reshape(alldeaths, [1,tot]));
    fprintf('\nmean birth rate = %5.3f, death rate = %5.3f.\n\n', [mbirth, mdeath])
    
    % variable wise summary
    P = size(betamat, 2);
    mat = zeros(P, 3);
    for i = 1:P
        betavec = betamat(:,i); betavec = betavec(betavec~=0);
        mat(i,:) = [mean(betavec), quantile(betavec,0.025), quantile(betavec, 0.975)];
    end
    a = mean(indmat); np = min(15, length(a));
    disp('posterior probability for important variables:')
    [a,I] = sort(a,'descend'); a = [I',mat(I,:), a']; disp(num2str(a(1:np,:))); fprintf('\n')
    % save('tmp.mat','a')
    
    if simu == 1
        disp('posterior probability for the true important variables:')
        disp(a(x0.inds)); fprintf('\n')
    end
    
    mat_var = a;
    
    gammas = reshape(allgammas, [1,tot]);
    plot(allgammas)
    disp('mean, lb and ub for gamma:')
    [l,u] = FindHPDset(gammas', 0.95, []); disp([mean(gammas),l,u]);
    
    tau2s = reshape(alltau2s, [1,tot]);
    plot(alltau2s)
    disp('mean, lb and ub for tau2:')
    [l,u] = FindHPDset(tau2s', 0.95, []); disp([mean(tau2s),l,u]);
    
    % model wise summary
    [xu, m, k] = unique(indmat, 'rows');
    count = histc(k, 1:size(xu, 1));
    [~,I] = sort(count,'descend');
    tab = [xu, count/sum(count)];
    tab = tab(I,:);
    % num2str(tab, 4)
    bmat = nan(size(tab,1),size(tab,2)-1);
    lb = bmat; ub = bmat;
    bmat0 = nan(size(tab,1),1); lb0 = bmat0; ub0 = bmat0;
    for i = 1:length(I)
        bmat(i,:) = mean(betamat(k==I(i), :),1);
        lb(i,:) = quantile(betamat(k==I(i), :),0.025,1);
        ub(i,:) = quantile(betamat(k==I(i), :),0.975,1);
        %[lb(i,:),ub(i,:)] = FindHPDset(betamat(k==I(i), :), 0.95, []);
        bmat0(i,:) = mean(beta0mat(k==I(i)));
        lb0(i,:) = quantile(beta0mat(k==I(i)),0.025);
        ub0(i,:) = quantile(beta0mat(k==I(i)),0.975);
    end
    
    save(strcat('paras',num2str(ncase),'.mat'),'allps', 'ps','tab','bmat','lb','ub','mat_var','bmat0','lb0','ub0')
end

disp('done')
end

function [LBout,UBout] = FindHPDset(Samples,p,npoints)
%Function to find the 100p % HPD set based on Samples
if isempty(npoints)
    npoints=200;
end

[f,x] = ksdensity(Samples,'npoints',npoints); N = size(Samples,1); maxf = max(f);
step = maxf/npoints;
HPDdensity = (step:step:maxf);
NHPD = size(HPDdensity,2);
LB = cell(1,NHPD); UB = cell(1,NHPD); Prob = zeros(1,NHPD);
for i=1:NHPD
    indices0 = find(HPDdensity(NHPD-i+1) < f);
    if ~isempty(indices0)
        indices1 = find(diff(indices0)> 1);
        if isempty(indices1)
            LB{i} = x(indices0(1)); UB{i} = x(indices0(end));
        elseif (size(indices1,1)==1)
            LB{i} = [x(indices0(1)) x(indices0(indices1(1)+1))];
            UB{i} = [x(indices0(indices1(1))) x(indices0(end))];
        else
            LB{i} = x(indices0(1)); UB{i} = [];
            for j=1:(size(indices1,2)-1)
                LB{i} = [LB{i} x(indices0(indices1(j)+1))];
                UB{i} = [UB{i} x(indices0(indices1(j)))];
            end
            UB{i} =[UB{i} x(indices0(end))];
        end
    end
    Ns = size(LB{i},2);
    count = zeros(1,Ns);
    for j=1:Ns
        count(j) = sum((LB{i}(j) <= Samples).*(Samples <= UB{i}(j)));
    end
    Prob(i) = sum(count)/N;
end
[minval indexmin] = min(abs(Prob - p));
LBout = LB{indexmin};
UBout = UB{indexmin};
end
