clear params;

VERSION = "m6-normalized-longer";

params.Y1 = 1/102;
params.Y2 = 1/102;
params.Y3 = 100/102;
%r11 = 0.005;
%r22 = 0.005;
%rs1 = 0.005;
%rs2 = 0.005;
r11 = (0.05:0.15:0.95);
r22 = (0.05:0.15:0.95);
rs1 = (0.05:0.15:0.95);
rs2 = (0.05:0.15:0.95);

%params.tspan = 200;
params.tspan = 10000;

paramsVector = [r11;r22;rs1;rs2];

for p=1:numel(r11)
    params.r11 = r11(p);
    for q=1:numel(r22)
        params.r22 = r22(q);
        for r=1:numel(rs1)
            params.rs1 = rs1(r);
            for s=1:numel(rs2)
                params.rs2 = rs2(s);
                data = modelSweepsNormalized(params);
                % record the densities
                alldata(p,q,r,s,:,:,:) = data.relab;
            end
        end
    end
end
save(strcat("data/",VERSION,"/relAbundance-",VERSION),"alldata")
save(strcat("data/",VERSION,"/params-", VERSION),"paramsVector")

%% Sweep across unique resource size of resource 1.
clear params;

VERSION = "m17";

Y1 = (0.1:0.5:10);
Y2 = 1 + zeros(1,20); 
Y3 = 1 + zeros(1,20);

paramsVector = [Y1;Y2;Y3];

for p=1:numel(Y1)    
    params.Y1 = Y1(p);
    for q = 1:numel(Y2)
        params.Y2 = Y2(q);
        for r=1:numel(Y3)
            params.Y3 = Y3(r);
            data = modelSweeps(params);
            % record the densities
            allData(p,q,r,:,:,:) = data.relab;
        end
    end
end
save(strcat("data/",VERSION,"/relAbundance-",VERSION),"allData")
save(strcat("data/",VERSION,"/params-", VERSION),"paramsVector")

%% Sweep across values of Y3 to recapitulate main text model.
clear params;

VERSION = "m47";

gamma = (0.025:0.025:1);
paramsVector = [gamma];

for q=1:numel(gamma)
    disp(gamma(q));
    params.Y1 = (1-gamma(q))/(2-gamma(q));
    params.Y2 = (1-gamma(q))/(2-gamma(q));
    params.Y3 = gamma(q)/(2-gamma(q));
    
    data = modelSweepsNormalized(params); 

    % record the densities
    alldata(q,:,:,:) = data.relab;
end
save(strcat("data/",VERSION,"/relAbundance-",VERSION),"alldata")
save(strcat("data/",VERSION,"/params-", VERSION),"paramsVector")