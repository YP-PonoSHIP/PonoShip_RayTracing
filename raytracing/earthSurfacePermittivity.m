function [epsilon, sigma, complexEpsilon] = earthSurfacePermittivity(mtl, fc, varargin)
%

% Copyright 2019-2024 The MathWorks, Inc.

%#codegen 

% Validate material name input
mtl = convertCharsToStrings(mtl);
validateattributes(mtl, {'string'}, {'scalar'}, ...
    'earthSurfacePermittivity', 'material name');

% Materials' frequency and temperature limits
mtlLib = {...
    % Name              Fmax    Tmin    Tmax      
    'pure-water'        1e12    -4      40;
    'sea-water'         1e12    -4      40;
    'pure-ice'          1e12    -60     0;
    'wet-ice'           1e12    NaN     NaN;
    'multi-year-ice'    1e11    -30     -2;
    'dry-snow'          1e11    -60     0;
    'wet-snow'          1e11    -60     0;
    'soil'              1e12    -273.15 Inf;
    'vegetation'        1e12    -20     Inf};
mtlNames = {mtlLib{:,1}};
mtlParams = cell2mat(reshape({mtlLib{:,2:4}}, [], 3));

% Use lower-case material name and set dry-ice = pure-ice
if strcmpi(mtl, "dry-ice")
    mtlLower = "pure-ice";
else
    mtlLower = lower(mtl);
end

% Validate material name
mtlMatchingIdx = find(strcmp(mtlLower, mtlNames));
coder.internal.errorIf(isempty(mtlMatchingIdx),...
    'shared_channel:earthSurfacePermittivity:InvalidMaterial', mtl);
mtlMatchingIdx = mtlMatchingIdx(1);

% Validate frequency
validateattributes(fc, {'numeric'}, ...
    {'real','positive','finite','scalar'}, ...
    'earthSurfacePermittivity', 'carrier frequency');
fc = double(fc);
Fmax = mtlParams(mtlMatchingIdx, 1);
coder.internal.errorIf(fc > Fmax, ...
    'shared_channel:earthSurfacePermittivity:FcOutOfRange', ...
    sprintf('%g', Fmax/1e9), mtl);
fcGHz = fc/1e9;

% Check number of inputs
coder.internal.errorIf( ...
    any(strcmp(mtlLower, {'pure-water','dry-ice','pure-ice','wet-ice'})) && (nargin ~=3), ...
    'shared_channel:earthSurfacePermittivity:IncorrectNumInputs', ...
    "3", mtl);
coder.internal.errorIf( ...
    any(strcmp(mtlLower, {'dry-snow','multi-year-ice','sea-water','vegetation'})) && (nargin ~= 4), ...
    'shared_channel:earthSurfacePermittivity:IncorrectNumInputs', ...
    "4", mtl);
coder.internal.errorIf( ...
    any(strcmp(mtlLower, {'wet-snow'})) && (nargin ~= 5), ...
    'shared_channel:earthSurfacePermittivity:IncorrectNumInputs', ...
    "5", mtl);
coder.internal.errorIf( ...
    strcmp(mtlLower, 'soil') && (nargin ~= 7) && (nargin ~= 8), ...
    'shared_channel:earthSurfacePermittivity:IncorrectNumInputs', ...
    "7 or 8", mtl);

% Validate temperature
T = coder.nullcopy(0);
Tmin = mtlParams(mtlMatchingIdx, 2);
Tmax = mtlParams(mtlMatchingIdx, 3);
if ~isnan(Tmin)
    validateattributes(varargin{1}, {'numeric'}, {'real','scalar'}, ...
        'earthSurfacePermittivity', 'temperature');
    T = double(varargin{1});

    if isfinite(Tmax)
        coder.internal.errorIf( T < Tmin || T > Tmax, ...
            'shared_channel:earthSurfacePermittivity:TempOutOfRange', ...
            sprintf('%g', Tmin), sprintf('%g', Tmax), mtl);
    else
        coder.internal.errorIf( T < Tmin, ...
            'shared_channel:earthSurfacePermittivity:TempLow', ...
            sprintf('%g', Tmin), mtl);
    end
end

% Calculate
complexEpsilon = coder.nullcopy(0);
switch mtlLower
    case 'pure-water'
        % Calculate permittivity
        complexEpsilon = getPureWaterPermittivity(fcGHz, T);
    case 'sea-water'
        if nargin == 4 % For codegen
            % Validate salinity in g/kg or ppt
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','finite','scalar','>=',0,'<=',40}, ...
                'earthSurfacePermittivity', 'sea water salinity');
            S = double(varargin{2});

            % Calculate permittivity. Sea water is equal to pure water when
            % salinity equals zero.
            if S == 0
                complexEpsilon = getPureWaterPermittivity(fcGHz, T);
            else
                complexEpsilon = getSeaWaterPermittivity(fcGHz, T, S);
            end
        end
    case 'pure-ice'
        % Calculate permittivity
        complexEpsilon = getPureIcePermittivity(fcGHz, T);
    case 'wet-ice' % T = 0
        % Validate liquid water volume fraction
        validateattributes(varargin{1}, {'numeric'}, ...
            {'real','scalar','>=',0,'<=',1}, ...
            'earthSurfacePermittivity', ...
            'wet ice liquid water volume fraction');
        F_wc = double(varargin{1});

        % Calculate permittivity
        complexEpsilon = getWetIcePermittivity(fcGHz, F_wc);
    case 'multi-year-ice'
        if nargin == 4 % For codegen
            % Validate air volume fraction
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',1}, ...
                'earthSurfacePermittivity', ...
                'multi-year ice air volume fraction');
            v_a = double(varargin{2});

            % Calculate permittivity. Multi-year ice is equal to pure ice when
            % the air volume fraction is zero.
            if v_a == 0
                complexEpsilon = getPureIcePermittivity(fcGHz, T);
            else
                complexEpsilon = getMultiYearIcePermittivity(fcGHz, T, v_a);
            end
        end
    case 'dry-snow'
        if nargin == 4 % For codegen
            % Validate dry snow density in g/cm^3
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','scalar','nonnegative','finite'}, ...
                'earthSurfacePermittivity', ...
                'dry snow density');
            rho_ds = double(varargin{2});

            % Calculate permittivity
            complexEpsilon = getDrySnowPermittivity(fcGHz, T, rho_ds);
        end
    case 'wet-snow'
        if nargin == 5 % For codegen
            % Validate dry snow density in g/cm^3
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','scalar','nonnegative','finite'}, ...
                'earthSurfacePermittivity', ...
                'wet snow density');
            rho_ds = double(varargin{2});

            % Validate liquid water volume fraction
            validateattributes(varargin{3}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',1}, ...
                'earthSurfacePermittivity', ...
                'snow liquid water volume fraction');
            F_wc = double(varargin{3});

            % Calculate permittivity
            complexEpsilon = getWetSnowPermittivity(fcGHz, T, rho_ds, F_wc);
        end
    case 'soil'
        if nargin >= 7 % For codegen
            % Validate sand percentage
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',100}, ...
                'earthSurfacePermittivity', 'sand percentage in soil');
            sandPercent = double(varargin{2});

            % Validate clay percentage
            validateattributes(varargin{3}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',100}, ...
                'earthSurfacePermittivity', 'clay percentage in soil');
            clayPercent = double(varargin{3});

            % Validate sand and clay percentage sum no larger than 1
            coder.internal.errorIf(clayPercent + sandPercent > 100, ...
                'shared_channel:earthSurfacePermittivity:InvalidSoilPercentage');

            % Validate specific gravity
            validateattributes(varargin{4}, {'numeric'}, ...
                {'real','finite','scalar','nonnegative'}, ...
                'earthSurfacePermittivity', 'soil specific gravity');
            rho_s = double(varargin{4});

            % Validate volumetric water content
            validateattributes(varargin{5}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',1}, ...
                'earthSurfacePermittivity', 'soil volumetric water content');
            m_v = double(varargin{5});

            % Get bulk density
            if nargin == 8
                % Validate bulk density input
                validateattributes(varargin{6}, {'numeric'}, ...
                    {'real','finite','scalar','nonnegative'}, ...
                    'earthSurfacePermittivity', 'soil bulk density');
                rho_b = double(varargin{6});
            else % Calculate bulk density using [ITU-R P.527-6, eq. (57)]
                % Ignore a constituent it its percentage if < 1%, but sum
                % up the included consitutent percentages to 100% [ITU-R
                % P.527-6, p. 18]
                siltPercent = (100 - sandPercent - clayPercent);
                P = [sandPercent; clayPercent; siltPercent];
                P(P<1) = 0;
                P = P.*(100/sum(P));
                P(P<1) = 1; % log(1) = 0

                mult = [0.078886; 0.038753; 0.032732];
                rho_b = 1.07256 + sum(mult.*log(P));
            end

            % Calculate permittivity
            complexEpsilon = getSoilPermittivity( ...
                fcGHz, T, sandPercent, clayPercent, rho_s, m_v, rho_b);
        end
    otherwise %  'vegetation'
        if nargin == 4 % For codegen
            % Validate gravimetric water content
            validateattributes(varargin{2}, {'numeric'}, ...
                {'real','scalar','>=',0,'<=',0.7}, ...
                'earthSurfacePermittivity', ...
                'vegetation gravimetric water content');
            M_g = double(varargin{2});

            % Calculate permittivity depending on temperature
            if T > 0
                complexEpsilon = ...
                    getVegetationPermittivityAboveFreezing(fcGHz, T, M_g);
            else
                complexEpsilon = ...
                    getVegetationPermittivityBelowFreezing(fcGHz, T, M_g);
            end
        end
end

% Real relative permittivity
epsilon = real(complexEpsilon);

% Absolute permittivity for free space (electric constant)
epsilon0 = 8.854187817e-12; % From [1]

% Derive conductivity from complex relative permittivity
sigma = -imag(complexEpsilon)*(2*pi*fc*epsilon0);

% Conductivity should be non-negative
coder.internal.errorIf(sigma < 0, 'shared_channel:earthSurfacePermittivity:InvalidConductivity')

end

function [eps_s, eps_1, eps_inf, f_1, f_2] = getWaterParams(T)
% ITU-R P.527-6, Section 5.1.1

theta = 300/(T + 273.15) - 1;               % Eq. (11)        
eps_s = 77.66 + 103.3*theta;                % Eq. (8) 
eps_1 = 0.0671*eps_s;                       % Eq. (9)
eps_inf = 3.52 - 7.52*theta;                % Eq. (10)
f_1 = 20.20 - 146.4*theta + 316*theta^2;    % Eq. (12)
f_2 = 39.8*f_1;                             % Eq. (13)

end

function [eps, eps_r, eps_c, f_1] = getPureWaterPermittivity(fcGHz, T)
% ITU-R P.527-6, Section 5.1.1

[eps_s, eps_1, eps_inf, f_1, f_2] = getWaterParams(T);

% Complex relative permittivity
eps_r = ...                             % Eq. (6)
    (eps_s - eps_1)   / (1 + (fcGHz/f_1)^2) + ...
    (eps_1 - eps_inf) / (1 + (fcGHz/f_2)^2) + ...
    eps_inf;                            
eps_c = ...                             % Eq. (7)
    (fcGHz/f_1) * (eps_s - eps_1)   / (1 + (fcGHz/f_1)^2) + ...
    (fcGHz/f_2) * (eps_1 - eps_inf) / (1 + (fcGHz/f_2)^2); 
eps = eps_r - 1i*eps_c;                 % Eq. (5)

end

function eps = getSeaWaterPermittivity(fcGHz, T, S)
% ITU-R P.527-6, Section 5.1.2

% Supplemental terms
[eps_s, eps_1, eps_inf, f_1, f_2] = getWaterParams(T);
eps_ss = eps_s * ...                    % Eq. (17)
    exp(-3.33330e-3*S + 4.74868e-6*S^2);
f_1s = f_1 * ...                        % Eq. (18) in GHz
    (1+S*(2.3232e-3 - 7.9208e-5*T + 3.6764e-6*T^2 + ...
    3.5594e-7*T^3 + 8.9795e-9*T^4));
eps_1s = eps_1 * ...                    % Eq. (19)
    exp(-6.28908e-3*S + 1.76032e-4*S^2 - 9.22144e-5*T*S);
f_2s = f_2 * ...                        % Eq. (20) in GHz
    (1 + S*(-1.99723e-2 + 1.81176e-4*T));
eps_infs = eps_inf * ...                % Eq. (21)
    (1 + S*(-2.04265e-3 + 1.57883e-4*T));    
sigma_35 = ...                          % Eq. (23)
    2.903602 + 8.607e-2*T + 4.738817e-4*T^2 - ...
    2.991e-6*T^3 + 4.3047e-9*T^4;
R_15 = S * ...                          % Eq. (24)
    (37.5109 + 5.45216*S + 1.4409e-2*S^2)/(1004.75 + 182.283*S + S^2);
alpha_0 = ...                           % Eq. (26)
    (6.9431 + 3.2841*S - 9.9486e-2*S^2) / (84.850 + 69.024*S + S^2);
alpha_1 = ...                           % Eq. (27)
    49.843 - 0.2276*S + 0.198e-2 * S^2; 
R_T15 = ...                             % Eq. (25)
    1 + alpha_0*(T - 15)/(alpha_1 + T);
sigma_sw = sigma_35 * R_15 * R_T15;     % Eq. (22) in S/m

% Complex relative permittivity
eps_r = ...                             % Eq. (15)
    (eps_ss - eps_1s)   / (1 + (fcGHz/f_1s)^2) + ...
    (eps_1s - eps_infs) / (1 + (fcGHz/f_2s)^2) + ...
    eps_infs;                            
eps_c =...                              % Eq. (16)
    (fcGHz/f_1s) * (eps_ss - eps_1s)   / (1 + (fcGHz/f_1s)^2) + ...
    (fcGHz/f_2s) * (eps_1s - eps_infs) / (1 + (fcGHz/f_2s)^2) + ...
    18 * sigma_sw / fcGHz;
eps = eps_r - 1i*eps_c;                 % Eq. (14)

end

function eps = getPureIcePermittivity(fcGHz, T)
% ITU-R P.527-6, Section 5.1.3.1

% Supplemental terms
theta = 300/(T + 273.15) - 1;           % Eq. (34)  
A = ...                                 % Eq. (31)
    (0.00504 + 0.0062*theta)*exp(-22.1*theta); 
tau = 335 / (T + 273.15);               % Eq. (33)
B =  ...                                % Eq. (32)
    0.0207 / (T + 273.15) * exp(-tau) / (exp(-tau)-1)^2 + ...
    1.16e-11 * fcGHz^2 + exp(-9.963 + 0.0372*T);

% Complex relative permittivity
eps_r = 3.1884 + 0.00091*T;             % Eq. (29)
eps_c = A / fcGHz + B * fcGHz;          % Eq. (30)
eps = eps_r - 1i*eps_c;                 % Eq. (28)

end

function eps = getWetIcePermittivity(fcGHz, F_wc)
% ITU-R P.527-5, Section 5.1.3.2; not in ITU-R P.527-6

eps_ice = getPureIcePermittivity(fcGHz, 0);
eps_pw  = getPureWaterPermittivity(fcGHz, 0);

% Complex relative permittivity
eps = ...                               % Eq. (35)
    (((eps_ice + 2*eps_pw) + 2*(eps_ice-eps_pw)*(1-F_wc)) / ...
     ((eps_ice + 2*eps_pw) -   (eps_ice-eps_pw)*(1-F_wc))) * eps_pw;
 
end

function eps = getMultiYearIcePermittivity(fcGHz, T, v_a)
% ITU-R P.527-6, Section 5.1.3.3.2

% Supplemental terms [p. 15]
eps_ice = getPureIcePermittivity(fcGHz, T);

% Complex relative permittivity [eq. (49)-(50)]
A = 2;
B = 1 - 2.*eps_ice - (3*v_a).*(1 - eps_ice);
C = -eps_ice;
eps = (-B + sqrt(B^2 - 4*A*C)) / (2*A); % Reference has wrong sign for sqrt

% Imaginary permittivity can exceed zero due to numerical precision
eps = real(eps) + 1i*min(imag(eps), 0);

end

function eps = getDrySnowPermittivity(fcGHz, T, rho_ds)
% ITU-R P.527-6, Section 5.1.4.1

% Supplemental terms [p. 15]
eps_ice = getPureIcePermittivity(fcGHz, T);
rho_ice = 0.916; % g/cm^3
f_ice = rho_ds / rho_ice;

% Complex relative permittivity [eq. (51)-(53)]
if rho_ds <= 0.5
    eps_r = 1 + 1.9*rho_ds;
else
    eps_r = 0.51 + 2.88*rho_ds;
end
eps_i = -3*imag(eps_ice)*f_ice * eps_r^2 * (2*eps_r+1) / ...
    ((real(eps_ice) + 2*eps_r) * (real(eps_ice) + 2*eps_r^2));
eps = eps_r - 1i*eps_i;

end

function eps = getWetSnowPermittivity(fcGHz, T, rho_ds, F_wc)
% ITU-R P.527-6, Section 5.1.4.2

% Supplemental terms [p. 16]
eps_ds = getDrySnowPermittivity(fcGHz, T, rho_ds);
[eps_pw, ~, ~, ~] = getPureWaterPermittivity(fcGHz, T);

% Complex relative permittivity [eq. (54)-(55)]
A = 2;
B = eps_pw - 2*eps_ds - 3*F_wc*(eps_pw - eps_ds);
C = -eps_pw*eps_ds;
eps = (-B + sqrt(B^2 - 4*A*C)) / (2*A);

% Imaginary permittivity can exceed zero due to numerical precision
eps = real(eps) + 1i*min(imag(eps), 0);

end

function eps = getSoilPermittivity(fcGHz, T, P_sand, P_clay, rho_s, m_v, rho_b)
% ITU-R P.527-6, Section 5.2

% Complex relative permittivity of free water [eq. (65)-(70)]
sigma_1 = 0.0467 + 0.2204*rho_b - 0.004111*P_sand - 0.006614*P_clay; 
sigma_2 = -1.645 + 1.939*rho_b - 0.0225622*P_sand + 0.01594*P_clay;
sigma_eff_r = (fcGHz/1.35) * (sigma_1 - sigma_2) / (1 + (fcGHz/1.35)^2);
sigma_eff_i = sigma_2 + (sigma_1 - sigma_2) / (1 + (fcGHz/1.35)^2);
[~, eps_pw_r, eps_pw_c] = getPureWaterPermittivity(fcGHz, T);
temp = 18 * (rho_s - rho_b) / (fcGHz * rho_s * m_v);
eps_fw_r = eps_pw_r + temp*sigma_eff_r;
eps_fw_i = eps_pw_c + temp*sigma_eff_i;

% Complex relative permittivity [eq. (58)-(64)]
alpha = 0.65;
beta_r = 1.2748 - 0.00519 * P_sand - 0.00152 * P_clay; 
beta_i = 1.33797 - 0.00603 * P_sand - 0.00166 * P_clay;
eps_sm_r = (1.01 + 0.44*rho_s)^2 - 0.062;
eps_r = (1 + (rho_b/rho_s)*(eps_sm_r^alpha - 1) + ...
    (m_v^beta_r) * (eps_fw_r^alpha) - m_v) ^ (1/alpha);
eps_i = ((m_v^beta_i) * (eps_fw_i^alpha)) ^ (1/alpha);
eps = eps_r - 1i*eps_i;

end

function eps = getVegetationPermittivityAboveFreezing(fcGHz, T, M_g)
% ITU-R P.527-6, Section 5.3.1

% Supplemental terms [eq. (75)-(77)]
eps_dv = 1.7 - 0.74*M_g + 6.16*M_g^2;
v_fw = M_g*(0.55*M_g - 0.076);
v_bw = 4.64*M_g^2/(1 + 7.36*M_g^2);

% Complex relative permittivity [eq. (72)-(74)]
[~, eps_pw_r, eps_pw_i, f_1] = getPureWaterPermittivity(fcGHz, T);
temp1 = sqrt(fcGHz/(0.02*f_1));
temp2 = fcGHz/(0.01*f_1);
eps_r = eps_dv + v_fw*eps_pw_r + ...
    v_bw*(2.9 + 55*(1 + temp1)/(1 + 2*temp1+ temp2));
eps_i = v_fw*(eps_pw_i + 22.86/fcGHz) + ...
    v_bw*(55*temp1/(1 + 2*temp1+ temp2));
eps = eps_r - 1i*eps_i;

end

function eps = getVegetationPermittivityBelowFreezing(fcGHz, T, M_g)
% ITU-R P.527-6, Section 5.3.2

% Supplemental terms [eq. (80)-(89)]
T_f = -6.5; % [p. 24]    
Delta = T - T_f;
c1 = fcGHz/1.2582;
c2 = 0.2054;
c2a = c2*pi/2;
c3 = 0.4108;
denom = 1 + 2*c1^c2 * cos(c2a) + c1^c3;
X1 = (1 + c1^c2 * cos(c2a)) / denom;
Y1 = c1^c2 * sin(c2a) / denom;
A_ice = 0.001 - 0.012*M_g + 0.0082*M_g^2;
B_ice = 0.036 - 0.2389*M_g + 0.1435*M_g^2; 
C_ice = -0.0538 + 0.4616*M_g - 0.3398*M_g^2;
v_ice = A_ice*Delta^2 + B_ice*Delta + C_ice;
eps_dv = 6.76 - 10.24*M_g + 6.19*M_g^2;
v_fw = (-0.106 + 0.6591*M_g - 0.610*M_g^2) * ...
    exp((0.06 + 0.6883*M_g + 0.0001*M_g^2)*Delta);
v_bw = (-0.16 + 1.1876*M_g - 0.387*M_g^2) * ...
    exp((0.721 - 1.2733*M_g + 0.8139*M_g^2)*Delta);

% Complex relative permittivity [eq. (72) & (78)-(79)]
c4 = fcGHz/9;
c5 = 82.2 / (1 + c4^2);
eps_r = eps_dv + v_fw*(4.9 + c5) + v_bw*(8.092 + 14.2067*X1) + 3.15*v_ice;
eps_i = v_fw*(c4*c5 + 11.394/fcGHz) + 14.2067*v_bw*Y1;
eps = eps_r - 1i*eps_i;

end

% [EOF]