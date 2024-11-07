function f = funcNs2019(x, D)
% Calculates heat transfer effectiveness (NTU method) for a plate-fin heat exchanger
% Reference:
% Zhang, L., Hu, T., Yang, Z. et al. 
% Elite and dynamic opposite learning enhanced sine cosine algorithm for 
% application to plate-fin heat exchangers design problem. 
% Neural Comput & Applic (2021). https://doi.org/10.1007/s00521-021-05963-2
% Coded by Tianyu Hu (Thierry)
% Modified on 21-02-2020

% Inputs:
% x - array of design variables [Lh, Lc, H, t, n, l, Nh]
% D - parameter related to the design

% Constants and operating conditions
mh   = 0.8962;    % Hot side mass flow rate (kg/s)
mc   = 0.8296;    % Cold side mass flow rate (kg/s)
Th   = 513;       % Hot side inlet temperature (K)
Tc   = 277;       % Cold side inlet temperature (K)
Ph   = 1e5;       % Hot side inlet pressure (Pa)
Pc   = 1e5;       % Cold side inlet pressure (Pa)
Cph  = 1017.7;    % Hot side specific heat (J/kg*K)
Cpc  = 1011.8;    % Cold side specific heat (J/kg*K)
rhoh = 0.8196;    % Hot side density (kg/m^3)
rhoc = 0.9385;    % Cold side density (kg/m^3)
muh  = 0.0000241; % Hot side dynamic viscosity (Pa*s)
muc  = 0.00002182; % Cold side dynamic viscosity (Pa*s)
Prh  = 0.6878;    % Hot side Prandtl number
Prc  = 0.6954;    % Cold side Prandtl number

% Extract design variables
Lh = x(1);     % Hot side length
Lc = x(2);     % Cold side length
H = x(3);      % Channel height
t = x(4);      % Fin thickness
n = x(5);      % Number of fins per unit length
l = x(6);      % Fin length
Nh = x(7);     % Number of hot-side fins

% Derived quantities
Hh = H;  % Hot side channel height
Hc = H;  % Cold side channel height
nh = n;  % Hot side fins per unit length
nc = n;  % Cold side fins per unit length
Ch = mh * Cph;  % Hot side capacity rate
Cc = mc * Cpc;  % Cold side capacity rate
Cmax = max(Ch, Cc);
Cmin = min(Ch, Cc);
Cr = Cmin / Cmax; % Heat capacity ratio

% Flow area and gap calculations
s = (1 / n - t); % Fin spacing
dh = (4 * s * (H - t) * l) / (2 * (s * l + (H - t) * l + t * (H - t)) + t * s); % Hydraulic diameter

% Cross-sectional areas
Affh = (Hh - t) * (1 - nh * t) * Lc * Nh;
Affc = (Hc - t) * (1 - nc * t) * Lh * (Nh + 1);

% Mass flux
Gh = mh / Affh;
Gc = mc / Affc;

% Reynolds numbers for hot and cold sides
Reh = (mh * dh) / (Affh * muh);
Rec = (mc * dh) / (Affc * muc);

% Colburn j-factors and friction factors for hot and cold sides
if Reh > 1500  % Turbulent flow (Re > 1500)
    jh = 0.21 * Reh^(-0.4) * (l / dh)^(-0.24) * (t / dh)^0.02;
    fh = 1.12 * Reh^(-0.36) * (l / dh)^(-0.65) * (t / dh)^0.17;
else            % Laminar flow (Re <= 1500)
    jh = 0.53 * Reh^(-0.5) * (l / dh)^(-0.15) * (s / H - t)^(-0.14);
    fh = 8.12 * Reh^(-0.74) * (l / dh)^(-0.41) * (s / H - t)^(-0.02);
end

if Rec > 1500
    jc = 0.21 * Rec^(-0.4) * (l / dh)^(-0.24) * (t / dh)^0.02;
    fc = 1.12 * Rec^(-0.36) * (l / dh)^(-0.65) * (t / dh)^0.17;
else
    jc = 0.53 * Rec^(-0.5) * (l / dh)^(-0.15) * (s / H - t)^(-0.14);
    fc = 8.12 * Rec^(-0.74) * (l / dh)^(-0.41) * (s / H - t)^(-0.02);
end

% Effective heat transfer areas
Ah = Lh * Lc * Nh * (1 + 2 * nh * (Hh - t));
Ac = Lh * Lc * (Nh + 1) * (1 + 2 * nc * (Hc - t));

% Pressure drop for the cold side
deltaPc = (2 * fc * Lc * (Gc^2)) / (rhoc * dh);

% NTU calculation
NTU = (1 / Cmin) * ((jh * Cph * Prh^(-2/3) * mh * Ah) / Affh + (jc * Cpc * Prc^(-2/3) * mc * Ac) / Affc);

% Heat exchanger effectiveness
epsilon = 1 - exp((1 / Cr) * NTU^0.22 * (exp(-Cr * NTU^0.78) - 1));

% Heat transfer calculation
Q = epsilon * Cmin * (Th - Tc);
Qx = abs(160000 - Q); % Difference from target heat transfer rate

% Total heat transfer area
AHT = Ah + Ac + 0.0001 * Qx;

% Output
f = AHT;  % Objective function value

end
