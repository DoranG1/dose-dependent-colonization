function data = model(params)

if isfield(params,'tspan')
    tspan = [0 params.tspan];
else
    tspan = [0 100];% From time 0 to 100
end

% Initial abundance of species (arbitrary example values)
if isfield(params,'X0')
    X0=params.X0;
else
    X0{1} = [0.001; 1; 1];
    X0{2} = [1; 0.001; 0.001];
end

% Initial resource vector
if isfield(params, 'Y1')
    Y1 = params.Y1;
else
    Y1 = 1;
end
if isfield(params, 'Y2')
    Y2 = params.Y2;
else
    Y2 = 1;
end
if isfield(params, 'Y3')
    Y3 = params.Y3;
else
    Y3 = 1;
end
if isfield(params, 'Y12')
    Y12 = params.Y12;
else 
    Y12 = 1;
end

if isfield(params, 'Y13')
    Y13 = params.Y13;
else
    Y13 = 1;
end

Y0 = [Y1 Y2 Y3 Y12 Y13]'; % Initial resource abundances.

if isfield(params,'r11')
    r11 = params.r11;
else
    r11 = 1;
end
if isfield(params,'r22')
    r22 = params.r22;
else
    r22 = 1;
end    
if isfield(params,'r33')
    r33 = params.r33;
else
    r33 = 1;    
end
if isfield(params,'rs1')
    rs1 = params.rs1;
else
    rs1 = 1;
end
if isfield(params,'rs2')
    rs2 = params.rs2;
else
    rs2 = 1;
end
if isfield(params,'rt1')
    rt1 = params.rt1;
else
    rt1 = 1;
end
if isfield(params,'rt3')
    rt3 = params.rt3;
else
    rt3 = 1;
end

R = [r11 0 0 rs1 rt1; % consumption rates of species 1
       0 r22 0 rs2 0; % consumption rates of species 2
       0 0 r33 0 rt3]; % consumption rates of species 3

% Concatenate initial conditions
if isfield(params,'numpass')
    numpass = params.numpass;
else
    numpass = 5;
    %numpass = 100;
end

if isfield(params,'dil')
    dil = params.dil;
else
    dil = 0.005;
end
if isfield(params,'makeplots')
    makeplots = params.makeplots;
else
    makeplots = 0;
end
for k=1:numel(X0)
    X = X0{k}*dil;
    for j=1:numpass
        Z0 = [X; Y0];
    
        % Use ode45 to simulate the system
        [T, Z] = ode45(@(t, Z) consumer_resource_model(t, Z, R), tspan, Z0);
        % display the relative abundance of each species
        relab = Z(end,1:3);
        %relab = Z(end,1:2);
        relab = relab/sum(relab);
        data.relab(k,j,:) = relab;

        disp([j relab])
        
        X = dil*Z(end,1:3)';
        %X = dil*Z(end,1:2)';
       
    end
end
end
% Define the ODE function
function dZdt = consumer_resource_model(t, Z, R)
    X = Z(1:3); % Extract species abundances
    %X = Z(1:2); % Extract species abundances
    Y = Z(4:end); % Extract resource levels
    
    dXdt = X .* (R * Y); % Change in species abundance
    dYdt = -Y .* (R' * X); % Change in resource levels
    
    dZdt = [dXdt; dYdt]; % Concatenate to form the derivatives vector
end