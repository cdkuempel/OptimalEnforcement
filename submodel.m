function TotalBiomass = submodel(p_PAs,Resolution)

% ============== Define model parameters ==============
NB = 0.7;     % Total enforcement budget
q = 0.7;      % Catchability 
d = 1;        % Market price of target species
cT = 0.018;   % Cost of travel of enforcers
cE = cT;      % Cost of travel of harvesters 
cB = 1;       % Penalty when caught
f = 2;        % Per-female fecundity
alpha = 1;    % Beverton-Holt parameter
beta = 0.01;  % Beverton-Holt parameter
z = 0.5;      % Diminishing returns exponent
m = 0.1;      % Natural mortality
M = 50;       % Total number of patches
P = p_PAs*M;  % Total number of protected patches


% Implement dilution of protection effort with increasing protected area extent
if p_PAs > 0 % If there are protected areas in the system, then find out how much each is enforced
   
   % This equation sets the probability of being caught as (NB - cT*P)./P
   % The enforcement budget per patch is the total budget, minus the cost of visiting each protected patch, divided by the total number of patchess:
   %        (NB - cT*P)/P
   % This is the probability of being caught, and is bounded between zero and one.
   
   pB = max(0,min(1,(NB-cT*P)./P));
   
else
   % If there are no protected areas, the chances of being caught poaching are zero (i.e., it's impossible).
   pB = 0;
end

% Generate all possible effort combinations
EffortVecP = linspace(0,4,Resolution);
EffortVecF = linspace(0,4,Resolution);
[Ep,Ef] = meshgrid(EffortVecP,EffortVecF); 

% The solution to the model is a double equilibrium - economic & ecological. 
% We need to find the effort combinations that yield such an equilibrium, 
% because that's where the system will end up.

% First work out what the populations (Nf, Np) would be if the system were at an economic equilibrium.
Nf = max(0,min(1,cT./(d.*q.*Ef.^(z-1))));
Np = max(0,min(1,(cT+cB.*pB)./(d.*q.*Ep.^(z-1))));

% At these population levels, the total larval arriving at each patch would be
L = f*((1-q.*Ef.^z).*(1-m).*Nf.*(M-P) + (1-q.*Ep.^z).*(1-m).*Np.*P)/M;

% Total number of settlers on any patch
S = alpha*L./(1+beta.*L); %Creates a matrix of number of settlers under different effort levels


% The following section runs one iteration of the population model for each 
% effort combination (inside and outside the reserves). We then look for situations
% where both timesteps have the same value (i.e., the ecological equilibrium).

Np_T = max(0,Np.*(1-q.*Ep.^z).*(1-m) + S); % Creates a matrix of biomass on protected patches in the next year based on varying harvesting effort levels
Nf_T = max(0,Nf.*(1-q.*Ef.^z).*(1-m) + S); % Creates a matrix of biomass on harvested (i.e. unprotected) patches in the next year based on varying harvesting effort levels

PopChange_P = Np_T - Np; % Creates a matrix of the change in biomass on protected patches under varying harvestinging effort levels
PopChange_F = Nf_T - Nf; % Creates a matrix of the change in biomass on harvested (i.e. unprotected) patches under varying harvestinging effort levels

Profit_P = d.*q.*Np_T.*Ep.^z - cT.*Ep - cB.*pB.*Ep; % Creates a matrix of the profit made on protected patches under varying harvesting effort levels
Profit_F = d.*q.*Nf_T.*Ef.^z - cT.*Ef; % Creates a matrix of the profit made on harvested (i.e. unprotected) patches under varying harvesting effort levels

Condition = abs(Profit_P) + abs(Profit_F) + abs(PopChange_F) + abs(PopChange_P);
% At very low levels of effort, the profits appear to be zero (since no one is fishing), but that's not the equilibrium we're looking for
Condition((Ef+Ep)<0.05) = inf; 

% The equilibrium point is one where profits are zero on both protected and unprotected reefs, and where the populations aren't changing.
[I,J] = find(Condition==min(Condition(:)));

TotalBiomass = P*Np_T(I,J) + (M-P)*Nf_T(I,J); %Total biomass should increase as number of protected areas and enforcement budget increase
