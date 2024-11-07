function f = funcNs2019(x, D)
% funcNs2019: NTU method for plate-fin heat exchanger optimization
%
% This function calculates the heat transfer effectiveness for a plate-fin heat exchanger.
% It employs the Number of Transfer Units (NTU) method to evaluate performance.
%
% Reference:
% Zhang, L., Hu, T., Yang, Z. et al. Elite and dynamic opposite learning enhanced sine cosine algorithm 
% for application to plate-fin heat exchangers design problem. Neural Comput & Applic (2021). 
% https://doi.org/10.1007/s00521-021-05963-2
%
% Coded by Tianyu Hu (Thierry)
% Modified on 21-02-2020
%
% Inputs:
%   x - Array of design variables [Lh, Lc, H, t, n, l, Nh]
%   D - Parameter related to the design

% Extract design variables
Lh = x(1); % Hot side length
Lc = x(2); % Cold side length
H  = x(3); % Channel height
t  = x(4); % Fin thickness
n  = x(5); % Number of fins per unit length
l  = x(6); % Fin length
Nh = x(7); % Number of hot-side fins

% Constants for operating conditions (Case 1)
mh   = 1.66;      % Hot side mass flow rate (kg/s)
mc   = 2;         % Cold side mass flow rate (kg/s)
Th   = 900 + 273.15; % Hot side inlet temperature (K)
Tc   = 200 + 273.15; % Cold side inlet temperature (K)
Ph   = 160000;    % Hot side inlet pressure (Pa)
Pc   = 200000;    % Cold side inlet pressure (Pa)
Cph  = 1122;      % Hot side specific heat (J/kg*K)
Cpc  = 1073;      % Cold side specific heat (J/kg*K)
rhoh = 0.6296;    % Hot side density (kg/m^3)
rhoc = 0.9638;    % Cold side density (kg/m^3)
muh  = 0.0000401; % Hot side dynamic viscosity (Pa*s)
muc  = 0.0000336; % Cold side dynamic viscosity (Pa*s)
Prh  = 0.731;     % Hot side Prandtl number
Prc  = 0.694;     % Cold side Prandtl number

% Heat capacity rates
Ch = mh * Cph;
Cc = mc * Cpc;
Cmax = max(Ch, Cc);
Cmin = min(Ch, Cc);
Cr = Cmin / Cmax; % Heat capacity ratio

% Derived quantities
Hh = H;   % Hot side channel height
Hc = H;   % Cold side channel height
nh = n;   % Hot side fins per unit length
nc = n;   % Cold side fins per unit length
s = (1 / n - t); % Fin spacing

% Hydraulic diameter based on 2019 formula
dh = (4 * s * (H - t) * l) / (2 * (s * l + (H - t) * l + t * (H - t)) + t * s);
Dhh = dh; % Hot side hydraulic diameter
Dhc = dh; % Cold side hydraulic diameter

Nc = Nh + 1; % Number of cold-side fins

% Cross-sectional areas for hot and cold sides
Affh = (Hh - t) * (1 - nh * t) * Lc * Nh; % Effective flow area (hot side)
Affc = (Hc - t) * (1 - nc * t) * Lh * Nc; % Effective flow area (cold side)

% Mass fluxes for hot and cold sides
Gh = mh / Affh;
Gc = mc / Affc;

% Reynolds numbers
Reh = (mh * dh) / (Affh * muh);
Rec = (mc * dh) / (Affc * muc);

% Colburn j-factors and friction factors for hot and cold sides
if Reh > 1500 % Turbulent flow
    jh = 0.21 * Reh^(-0.4) * (l / dh)^(-0.24) * (t / dh)^0.02;
    fh = 1.12 * Reh^(-0.36) * (l / dh)^(-0.65) * (t / dh)^0.17;
else          % Laminar flow
    jh = 0.53 * Reh^(-0.5) * (l / dh)^(-0.15) * (s / H - t)^(-0.14);
    fh = 8.12 * Reh^(-0.74) * (l / dh)^(-0.41) * (s / H - t)^(-0.02);
end

if Rec > 1500 % Turbulent flow
    jc = 0.21 * Rec^(-0.4) * (l / dh)^(-0.24) * (t / dh)^0.02;
    fc = 1.12 * Rec^(-0.36) * (l / dh)^(-0.65) * (t / dh)^0.17;
else          % Laminar flow
    jc = 0.53 * Rec^(-0.5) * (l / dh)^(-0.15) * (s / H - t)^(-0.14);
    fc = 8.12 * Rec^(-0.74) * (l / dh)^(-0.41) * (s / H - t)^(-0.02);
end

% Pressure drop calculations
deltaPh = (2 * fh * Lh * (Gh^2)) / (rhoh * Dhh);
deltaPc = (2 * fc * Lc * (Gc^2)) / (rhoc * Dhc);

% Heat transfer areas
Ah = Lh * Lc * Nh * (1 + 2 * nh * (Hh - t));
Ac = Lh * Lc * Nc * (1 + 2 * nc * (Hc - t));

% NTU calculation
NTU = (1 / Cmin) * ((jh * Cph * Prh^(-2/3) * mh * Ah) / Affh + (jc * Cpc * Prc^(-2/3) * mc * Ac) / Affc);

% Heat exchanger effectiveness
epsilon = 1 - exp((1 / Cr) * NTU^0.22 * (exp(-Cr * NTU^0.78) - 1));

% Heat transfer and objective function
Q = epsilon * Cmin * (Th - Tc);
Qx = abs(1070000 - Q); % Difference from target heat transfer rate
Ns = ((1 - epsilon) * ((Tc - Th)^2 / (Tc * Th)) + (deltaPh / Ph) + (deltaPc / Pc)) + (0.0001 * Qx);

% Output
f = Ns; % Objective function value

end
