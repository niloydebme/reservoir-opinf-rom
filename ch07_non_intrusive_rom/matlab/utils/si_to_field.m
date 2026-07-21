function fu = si_to_field()
%SI_TO_FIELD  Unit conversion factors: SI (MRST internal) → PTA/RTA field units
%
%  Usage:  fu = si_to_field();
%          p_psia = p_Pa * fu.Pa_to_psia;
%
%  All factors multiply the SI quantity to get the field-unit quantity.

% --- Pressure ---
fu.Pa_to_psia    = 1 / 6894.757;   % Pa      → psia
fu.Pa_to_psi     = 1 / 6894.757;   % Pa      → psi  (delta-p, same factor)

% --- Compressibility ---
fu.Pam1_to_psi1  = 6894.757;       % Pa⁻¹    → psi⁻¹

% --- Volumetric rates ---
fu.m3s_to_STBd   = 86400 / 0.158987; % m³/s  → STB/day  (oil or water)
fu.m3s_to_MSCFd  = 86400 / 28.3168;  % m³/s  → MSCF/day (gas)

% --- Permeability ---
fu.m2_to_mD      = 1 / 9.869233e-16; % m²    → mD

% --- Length / radius ---
fu.m_to_ft       = 3.28084;          % m      → ft

% --- Area ---
fu.m2_to_ft2     = 10.7639;          % m²     → ft²

% --- Viscosity ---
fu.Pas_to_cp     = 1e3;              % Pa·s   → cP

% --- Density ---
fu.kgm3_to_lbft3 = 0.0624279;       % kg/m³  → lb/ft³

% --- Time ---
fu.s_to_hr       = 1 / 3600;        % s      → hours
fu.day_to_hr     = 24;              % days   → hours

% --- Formation volume factor (dimensionless ratio — no conversion) ---
% Bo [RB/STB], Bg [rcf/scf]: user supplies directly, no SI equivalent in MRST

% --- Permeability-thickness ---
% kh [mD-ft] = k[mD] * h[ft]  (computed in Python from k_mD and h_ft)

fprintf('[si_to_field] Conversion factors loaded (SI → PTA/RTA field units)\n');
end