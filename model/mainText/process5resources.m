%% analyze what happens when you change the size of the unique niche
% define gamma as the fraction of the niche of species 1
% that overlaps with (i.e. is a shared resource with) species 2
% vary gamma from 2.5% to 100% and calculate alpha accordingly
gamma = ([0.025:0.025:1]);
%gamma(end) = 0.999;
%alpha = (1-gamma)./(2*gamma);
clear params;
params.dil = 1/200;
for q=1:numel(gamma)
    disp(gamma(q));
    params.gamma = gamma(q);
    params.X0{1} = [0.001;1;0];
    params.X0{2} = [0.01;1;0];
    params.X0{3} = [0.1;1;0];
    params.X0{4} = [1;1;0];
    params.X0{5} = [1;0.1;0];
    params.X0{6} = [1;0.01;0];
    params.X0{7} = [1;0.001;0];
    
    data = model5resources(params); 

    % record the densities
    alldata(q,:,:,:) = data.relab;
end

%save("relAbundanceByGamma","alldata")
save("relAbundanceByGamma-5resources","alldata")

%% vary r21 and calculate relative abundances at all ratios

clear params;

gamma = ([0.025:0.025:1]);
%alpha = (1-gamma)./(2*gamma);
r21 = [1/1.2, 1/1.15, 1/1.1];
params.dil = 1/200;
params.numpass = 50;

for r=1:numel(r21)
    for x=1:numel(gamma)
        params.gamma = gamma(x);
        params.r21 = r21(r);
        params.X0{1} = [0.001;1;1];
        params.X0{2} = [0.01;1;1];
        params.X0{3} = [0.1;1;1];
        params.X0{4} = [1;1;1];
        params.X0{5} = [1;0.1;1];
        params.X0{6} = [1;0.01;1];
        params.X0{7} = [1;0.001;1];
        
        data = model5resources(params); 

        % record the densities
        alldata(x,r,:,:,:) = data.relab;
    end
end

save("relAbundanceByR21GammaAsymmetric-5resources","alldata")

%% let's test the betas across different values of gamma
% define gamma13 as the fraction of the niche of species 1
% that is shared with species 3
% gamma13 = 2*alpha*beta / (2*alpha*(1-beta)+1)
% beta = gamma13 * (2*alpha+1)/(2*alpha*(1+gamma13))
% gamma13 gives some weird behaviors, so stick to defining beta
gamma = ([0.025:0.025:1]);
%gamma(end) = 0.999;
%alpha = (1-gamma)./(2*gamma);
beta = [0 0.05 0.1 0.25 0.5 0.75 0.8 0.85 0.9 0.975];
clear params;

for q=1:numel(beta)
    for x=1:numel(gamma)
        params.gamma=gamma(x);
        disp(beta(q));
        params.r3 = 10;
        %params.r3 = 1;
        %params.c = 1;
        params.c = (0.1)/(1.1); % set c to unique resource size when gamma=0.9
        params.beta = beta(q);
        % Species 3 starting abundance is equal to that of species 2.
        params.X0{1} = [0.001;1;1];
        params.X0{2} = [1;0.001;0.001];
        %data = model(params);
        data = model5resources(params);
    
        % record the densities
        alldataBetaGamma(x, q,:,:,:) = data.relab;
    end
end

%save("relAbundanceByBetaGammaAll","alldataBetaGamma")
save("relAbundanceByBetaGammaAll-5resources","alldataBetaGamma")

%% Test how the relative consumption of the shared resource by species 1 and 2
% affects the degree of dose dependence
% Keep r22 (the consumption rate of the shared resource by species 2) at 1.
% Vary r21 (the consumption rate of the shared resource by species 1).
clear params;
gamma = ([0.025:0.025:1]);
%alpha = (1-gamma)./(2*gamma);
r21 = [0.9, 1, 1.1];
% Test a range of r21 with more skewed dose dependences.
%r21 = [0.1, 0.25, 0.5];

for r=1:numel(r21)
    for x=1:numel(gamma)
        params.gamma=gamma(x);
        params.r21 = r21(r);
    
        % record the densities
        data = model5resources(params);
        alldataR21Gamma(x, r,:,:,:) = data.relab;
    end
end

save("relAbundanceByR21GammaAll-5resources","alldataR21Gamma")
%save("relAbundanceByR21GammaAsymmetric","alldataR21Gamma")

%% Test how the relative consumption of the shared resource by species 1 and 3
% affects the degree of dose dependence
% Vary r3 (the consumption rate of the shared resource by species 3).
gamma = ([0.025:0.025:1]);
%alpha = (1-gamma)./(2*gamma);
r3 = 10.^([-2:0.5:2]);
beta = 0.5;
clear params;
for r=1:numel(r3)
    for x=1:numel(gamma)
        params.gamma=gamma(x);
        params.r3 = r3(r);
        params.r13 = 1;
        params.beta = beta;
        params.c = (0.1)/(1.1); % set c to unique resource size when gamma=0.9
        % Species 3 starting abundance is equal to that of species 2.
        params.X0{1} = [0.001;1;1];
        params.X0{2} = [1;0.001;0.001];
    
        % record the densities
        data = model5resources(params);
        alldataR3Gamma(x, r,:,:,:) = data.relab;
    end
end

save("relAbundanceByR3GammaAll-5resources","alldataR3Gamma")

%% analyze what happens when you change the size of the unique niche
% define gamma as the fraction of the niche of species 1
% that overlaps with (i.e. is a shared resource with) species 2
% vary gamma from 2.5% to 100% and calculate alpha accordingly
% let the model run to 100 passages to allow the communities to equlibrate.
gamma = ([0.025:0.025:1]);
%gamma(end) = 0.999;
clear params;
params.dil = 1/200;
params.numpass = 75;
for q=1:numel(gamma)
    disp(gamma(q));
    params.gamma = gamma(q);
    params.X0{1} = [0.001;1;0];
    params.X0{2} = [1;0.001;0];
    
    data = model5resources(params); 

    % record the densities
    alldataMorePassages(q,:,:,:) = data.relab;
end

save("relAbundanceByGammaMorePassages-5resources","alldataMorePassages")


