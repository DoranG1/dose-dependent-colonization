function data = model5resources(params)

if isfield(params,'tspan')
    tspan = [0 params.tspan];
else
    %tspan = [0 100];% From time 0 to 100
    tspan = [0 10000];
end

% Initial abundance of species (arbitrary example values)
%X0 = [1; 1; 1]; % 3 species
%X0 = [1; 0.001; 1];
if isfield(params,'X0')
    X0=params.X0;
else
    X0{1} = [0.001; 1; 0];
    X0{2} = [1; 0.001; 0];
end

% Initial abundance of resources (arbitrary example values)
%if isfield(params,'alpha')
%    alpha = params.alpha;
%else
%    alpha = 0.05;
%end

if isfield(params,'gamma')
    gamma = params.gamma;
else
    gamma = 0.05;
end

% ratio of competitor in species 3 to species 1 niche
if isfield(params,'beta')
    beta = params.beta;
else
    % CHANGELOG: beta is now 0 after flipping the parameters in Y0 below.
    beta = 0;
end

if isfield(params, 'c')
    c = params.c;
else
    c = 0;
end

% CHANGELOG: SWAPPED 1-beta and beta
%Y0 = [2*alpha*(1-beta) alpha 0. 1 2*alpha*beta alpha]'; % 6 resources
%Y0 = [2*alpha*(1-beta) alpha 2*alpha 1 2*alpha*beta alpha]'; % 6 resources
%Y0 = [(1-beta) 1 c 1/(2*alpha) beta]'; % 5 resources
%Y0 = [(1-beta) 1 c gamma/(1-gamma) beta]'; % un-normalized
%Y0 = [((1-gamma)*(1-beta))/(2-gamma) (1-gamma)/(2-gamma) (c*(1-gamma))/(2-gamma) gamma/(2-gamma) (beta*(1-gamma))/(2-gamma)]'; % normalized to the sum of a, b, ab
%Y0 = [((1-gamma)*(1-beta))/(2+c-gamma*(1+c)) (1-gamma)/(2+c-gamma*(1+c)) (c*(1-gamma))/(2+c-gamma*(1+c)) gamma/(2+c-gamma*(1+c)) (beta*(1-gamma))/(2+c-gamma*(1+c))]'; % normalized to the sum of all resources
Y0 = [((1-gamma)*(1-beta))/(2-gamma) (1-gamma)/(2-gamma) c gamma/(2-gamma) (beta*(1-gamma))/(2-gamma)]'; % normalized to the sum a, b, ab, fixed value of c

% Consumption rates matrix R (arbitrary example values)
%R = [1 0 0 1 1 0;
%    0 1 0 1 0 0;
%    0 0 1 0 1.5 0];
if isfield(params,'r11')
    r11 = params.r11;
else
    r11 = 1; %2;
end
if isfield(params,'r21')
    r21 = params.r21;
else
    r21 = 1;
end
if isfield(params,'r22')
    r22 = params.r22;
else
    r22 = 1;
end
if isfield(params,'r3')
    r3 = params.r3;
else
    r3 = 1;
end
if isfield(params,'r13')
    r13 = params.r13;
else
    r13 = r11;
end

R = [r11 0 0 r21 r13; % consumption rates of species 1
       0 1 0 r22 0; % consumption rates of species 2
       0 0 1 0   r3]; % consumption rates of species 3
%R = [0.1, 0.2, 0.3, 0.4, 0.1, 0.1;
%     0.2, 0.1, 0.4, 0.1, 0.2, 0.2;
%     0.3, 0.3, 0.1, 0.2, 0.1, 0.3];

% Concatenate initial conditions
if isfield(params,'numpass')
    numpass = params.numpass;
else
    numpass = 5;
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
        relab = relab/sum(relab);
        data.relab(k,j,:) = relab;

        disp([j relab])
        
        X = dil*Z(end,1:3)';
    
        % Plot the results
        if makeplots
            figure(1);
            subplot(2,1,1);
            plot(T, Z(:,1:3));
            title('Species Abundance');
            legend('Species 1', 'Species 2', 'Species 3');
        
            subplot(2,1,2);
            plot(T, Z(:,4:end));
            title('Resource Levels');
            legend('Resource 1', 'Resource 2', 'Resource 3', 'Resource 4', 'Resource 5', 'Resource 6');
            %pause;
        end
    end
    if makeplots
        % display the relative abundance of each species
        relab = Z(end,1:3);
        relab = relab/sum(relab);
        disp(relab)
    end
end
end
% Define the ODE function
function dZdt = consumer_resource_model(t, Z, R)
    X = Z(1:3); % Extract species abundances
    Y = Z(4:end); % Extract resource levels
    
    dXdt = X .* (R * Y); % Change in species abundance
    dYdt = -Y .* (R' * X); % Change in resource levels
    
    dZdt = [dXdt; dYdt]; % Concatenate to form the derivatives vector
end