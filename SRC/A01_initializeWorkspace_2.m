% Codi McKee
% Texas A&M University
% First Created: 18-FJan-2020
% Last Modified: 18-Jan-2020

% Initialize Workspace

if exist('Frames','var') == 0
    loadTimetemp = tic;
    fprintf('\n\n Loading Frames.mat...\n');
    load Frames.mat
    timeLoad = toc(loadTimetemp);
    fprintf('\n Loading COMPLETE!!\n');
    
elseif exist('Frames','var') == 1
    clearvars -except Frames Response FramesOUT_Pushover
else
    error('Variable "Frames" not in the workspace!')
end

addpath([pwd '\Input'])
addpath([pwd '\Material Models']);
addpath([pwd '\Output']);
addpath([pwd '\Other Functions']);
addpath([pwd '\Data\FEMA_P695']);
oldDir = cd;
cd([pwd '\RunFiles']);
newDir = cd;
cd(oldDir);

load Rebar.mat Rebar
load GroundMotions.mat
u = setUnits('SI');

tog.Buckling        = 1;
tog.Steel02         = 0;

if tog.Buckling
    nameSuffix = '';
elseif tog.Steel02
    nameSuffix = '_Steel02';
else
    nameSuffix = '_NoBuckling';
end

tRange          = [0 25 50 75 100];
frameNumRange   = [1];
numRange        = [34];
scaleRange      = [1];
ageVal          = 1;
corrosion_level = 2;


[Cs,Dc,Ccr,k,v_rs] = corrosionParams(corrosion_level,u);


%**************************************************************************
% Internal Functions
%

function [Cs,Dc,Ccr,k,v_rs] = corrosionParams(corrosion_level,u)
% Corrosion Parameters

if corrosion_level == 0
    Cs  = 0    * (u.kg/u.m3);     % [kg/m^3] Surface chloride concentration (marine environment)
elseif corrosion_level == 1
    Cs  = 1.75 * (u.kg/u.m3);     % [kg/m^3] Surface chloride concentration (marine environment)
elseif corrosion_level == 2
    Cs  = 2.95 * (u.kg/u.m3);     % [kg/m^3] Surface chloride concentration (marine environment)
end

Dc  = 1.29 * (u.cm2/u.yr);    % [cm/yr]  Diffusion coefficient
Ccr = 0.90 * (u.kg/u.m3);     % [kg/m^3] Critical chloride concentration
k   = 0.1;                    % Bar roughness factor
v_rs = 2;                     % Volumentric expansion ratio
% t   = 0  .* (u.yr);          % [yr] Time

end



