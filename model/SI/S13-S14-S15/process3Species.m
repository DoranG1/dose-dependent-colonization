%% sweep rt1 and rt3. 
clear params;

VERSION = "s1-m22";

params.Y1 = 1;
params.Y2 = 1;
params.Y3 = 10;
params.Y12 = 5;
params.Y13 = 40;
params.r11 = 1;
params.r22 = 1;
params.r33 = 1;
params.rs1 = 1;
params.rs2 = 1;
rt1 = (0.1:0.2:1);
rt3 = (0.1:0.2:1);

paramsVector = [rt1;rt3];

for p=1:numel(rt1)
    params.rt1 = rt1(p);
    for q=1:numel(rt3)
        params.rt3 = rt3(q);
        data = model3species(params);
        % record the densities
        alldata(p,q,:,:,:) = data.relab;
    end
end
save(strcat("data/",VERSION,"/relAbundance-",VERSION),"alldata")
save(strcat("data/",VERSION,"/params-", VERSION),"paramsVector")

%% sweep rs1, rs2, rt1, rt3.
clear params;

VERSION = "s2-m42-normalized-longer";

params.Y1 = 1/202;
params.Y2 = 1/202;
params.Y3 = 1/202;
params.Y12 = 100/202;
params.Y13 = 100/202;
r11 = 1;
r22 = 1;
r33 = 1;
rs1 = (0.05:0.15:0.95);
rs2 = (0.05:0.15:0.95);
rt1 = (0.05:0.15:0.95);
rt3 = (0.05:0.15:0.95);

params.tspan = 10000;

paramsVector = [rs1;rs2;rt1;rt3];

for p=1:numel(rs1)
    params.rs1 = rs1(p);
    for q=1:numel(rs2)
        params.rs2 = rs2(q);
        for r=1:numel(rt1)
            params.rt1 = rt1(r);
            for s=1:numel(rt3)
                params.rt3 = rt3(s);
                data = model3species(params);
                % record the densities
                alldata(p,q,r,s,:,:,:) = data.relab;
            end
        end
    end
end
save(strcat("data/",VERSION,"/relAbundance-",VERSION),"alldata")
save(strcat("data/",VERSION,"/params-", VERSION),"paramsVector")