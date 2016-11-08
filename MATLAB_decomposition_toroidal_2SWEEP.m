% Multipole Decomposotion with toroidal moment separation.
% ver 4.0

clc
clear all;

norm_length = 1e9;  % value to multiply for m -> nm converting  
                    % величина, на которую домножаем, чтобы перевести метры в нанометры
% norm_length = 1e6;    % value to multiply for m -> um converting  
%                       % величина, на которую домножаем, чтобы перевести метры в микрометры

n_min = 1;  % minimum graphic slider parameter and first considered parameter                   
            % минимальное значение слайдера; случаи с меньшим номером параметра не будут выводиться
n = 1;      % The original (first after fig creating) graphic slider parameter                  
            % первоначальное значение слайдера

FontSize = 16;          % font size at the titles % размер шрифта в подписях
FontSizeLeg = 9;       % font size at the legends % размер шрифта в легендах
LineWidth = 2.3;        % Line Width at the graphics % толщина линии на графиках            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent Variables
rad = 100e-9;   % base edge of pyramid/parallelepiped/etc or radius of sphere/cilinder/cone/etc, in METERS, MANUALLY, for Normalization! MUST BE CHECKED  
                % сторона основания пирамиды В МЕТРАХ, ПРОВЕРЯТЬ если нормируем на эффективное сечение!
geomCS = rad^2; % effective geometrical cross-section for normalization, if it has square shape 
                % эффективное поперечное сечение рассеяния для нормировки на него
%geomCS = pi*rad^2; % effective geometrical cross-section for normalization, if it has circle shape. 
% geomCS = 1e-18; % for usual normalization by nm  
%                 % нормировка на нм обычная!
% geomCS = 1e-12; % for usual normalization by um  
%                 % нормировка на мкм обычная!
E0x = 1 ; 
E0 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
mu0 = 1.25663706 * 10^(-6); % Vacuum permeability (magnetic constant)     % магнитная постоянная 
eps0 = 8.85 * 10^(-12);     % Vacuum permittivity (electric constant)     % электрическая постоянная
c = 2.998e+8;               % speed of light                              % скорость света
% Z = 120*pi;               % for phase diagrams, still legacy parameter  % для фазовых диаграмм
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fre              = dlmread ('fre.dat');
H                = dlmread ('H.dat');
Px               = dlmread ('Px.dat');
Py               = dlmread ('Py.dat');
Pz               = dlmread ('Pz.dat');
Tx               = dlmread ('Tx.dat');
Ty               = dlmread ('Ty.dat');
Tz               = dlmread ('Tz.dat');
mx               = dlmread ('mx.dat');
my               = dlmread ('my.dat');
mz               = dlmread ('mz.dat');
Mxx              = dlmread ('Mxx.dat');
Mxy              = dlmread ('Mxy.dat');
Mxz              = dlmread ('Mxz.dat');
Myx              = dlmread ('Myx.dat');
Myy              = dlmread ('Myy.dat');
Myz              = dlmread ('Myz.dat');
Mzx              = dlmread ('Mzx.dat');
Mzy              = dlmread ('Mzy.dat');
Mzz              = dlmread ('Mzz.dat');
Qxx              = dlmread ('Qxx.dat');
Qyy              = dlmread ('Qyy.dat');
Qzz              = dlmread ('Qzz.dat');
Qxy              = dlmread ('Qxy.dat');
Qxz              = dlmread ('Qxz.dat');
Qyx              = dlmread ('Qyx.dat');
Qyz              = dlmread ('Qyz.dat');
Qzx              = dlmread ('Qzx.dat');
Qzy              = dlmread ('Qzy.dat');
Oxxx             = dlmread ('Oxxx.dat');
Oxxy             = dlmread ('Oxxy.dat');
Oxxz             = dlmread ('Oxxz.dat');
Oyyx             = dlmread ('Oyyx.dat');
Oyyy             = dlmread ('Oyyy.dat');
Oyyz             = dlmread ('Oyyz.dat');
Ozzx             = dlmread ('Ozzx.dat');
Ozzy             = dlmread ('Ozzy.dat');
Ozzz             = dlmread ('Ozzz.dat');
Oxyz             = dlmread ('Oxyz.dat');
Lambdax          = dlmread ('Lambdax.dat');
Lambday          = dlmread ('Lambday.dat');
Lambdaz          = dlmread ('Lambdaz.dat');
absCS            = dlmread ('absCS.dat');
scat             = dlmread ('scat.dat');
ForScat          = dlmread ('ForScat.dat');
BackScat         = dlmread ('BackScat.dat');
ForScatPoint     = dlmread ('ForScatPoint.dat');
BackScatPoint    = dlmread ('BackScatPoint.dat');
ForScatPow       = dlmread ('ForScatPow.dat');
BackScatPow      = dlmread ('BackScatPow.dat');
ForScatPointPow  = dlmread ('ForScatPointPow.dat');
BackScatPointPow = dlmread ('BackScatPointPow.dat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_max = size(Px, 1);
epsd = 1;   % Permittivity of environment outside of the particle 
            % диэлектричекая проницаемость среды снаружи частицы
lambda = c./fre;    % wavelength as speed of light divided by frequency 
                    % длина волны - скорость света делить на частоту
k0 = - 2*pi*fre/c;  % k-vector % волновой вектор
                    % minus here, cause at COMSOL there is minus at the imaginary exponent e^(-i k*x), describing incident radiation
                    % минус потому что в COMSOL задан минус в мнимой экспоненте e^(-i k*x), описывающей падающее излучение

kd = k0 .* sqrt(epsd);  % k-vector in the environment outside of the particle 
                        % волновой вектор в среде
vd2 = c ./ sqrt(epsd); % speed of light in the environment outside of the particle, vd for scattering cross-section without dependence of minus in k0 
                       % скорость света в среде, vd для сечения рассеяния без зависимости от минуса в k0

% vd = c ./ sqrt(epsd);
vd = 2.*pi.*fre./( k0 .* sqrt(epsd) );  % speed of light in the environment too, with this vd extinction calculations works well
                                        % тоже скорость света в среде. С этим vd хорошо считается экстинкция.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extinction components for the wave, directed along Z and polarized along X
% Компоненты экстинкции для волны направленной по Z и поляризованной по X
ExtPx   = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* Px);
ExtTx   = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* ((1i.*kd./vd).*1.*Tx) );
ExtQxz  = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (-(1i.*kd./6).*Qxz) );
Extmy   = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (my./vd) );
ExtMyz  = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* ((1i.*kd./(2.*vd)).*Myz) );
ExtOxzz = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (-(kd.*kd./6) .* (Ozzx-Lambdax)) );


TxK = ((1i.*kd./vd) .* Tx);
TyK = ((1i.*kd./vd) .* Ty);
TzK = ((1i.*kd./vd) .* Tz);

% Extinction as sum of multipole components
% Экстинкция как сумма мультипольных составляющиx
ExtCS = ExtPx+ExtTx+ExtQxz+Extmy+ExtMyz+ExtOxzz;

% Components of Total Electric Dipole Moment (TED)
% Компоненты полного электрического дипольного момента
Dx = Px - (1i.*k0.*epsd./c).*Tx;
Dy = Py - (1i.*k0.*epsd./c).*Ty;
Dz = Pz - (1i.*k0.*epsd./c).*Tz;

% I incident - Poynting vector of incident wave, power flux density
% вектор пойнтинга падающей волны, плотность потока мощности
Iinc = ((eps0*epsd/mu0)^0.5*E0/2);

% Contributions of separate miltipoles at the power scattering, normalized by Poynting vector
% Вклады отдельных мультиполей в рассеяние мощности нормированные на вектор пойнтинга
ScatD = (k0.^4./(12.*pi.*eps0.^2.*vd2.*mu0)) .* (abs(Dx).^2 + abs(Dy).^2 + abs(Dz).^2) ./ Iinc;     % Electric Dipole % диполь

Scatm = ((k0.^4.*epsd./(12.*pi*eps0.*vd2)) .* (abs(mx).^2+abs(my).^2+abs(mz).^2)) ./ Iinc;  % Magnetic Dipole % магнитный диполь

ScatQ=((k0.^6.*epsd./(1440.*pi.*eps0.^2.*vd2.*mu0))...                                  % Electric Quadrupole % электрический квадруполь
    .*(abs(Qxx).^2+abs(Qxy).^2+abs(Qxz).^2+abs(Qyx).^2+abs(Qyy).^2+abs(Qyz).^2+...
    abs(Qzx).^2+abs(Qzy).^2+abs(Qzz).^2)) ./ Iinc;                                          

ScatM = ((k0.^6.*epsd.^2./(160.*pi*eps0.^1.*vd2)) ...                                   % Magnetic quadrupole % магнитный квадруполь
    .*(abs(Mxx).^2+abs(Mxy).^2+abs(Mxz).^2+abs(Myx).^2+abs(Myy).^2+abs(Myz).^2+...
    abs(Mzx).^2+abs(Mzy).^2+abs(Mzz).^2) ) ./ Iinc;                                         

% octupole parts from COMSOL % запчасти от октуполя из комсола
FullOxyz = Oxyz;
FullOxxy = Oxxy-Lambday;
FullOxxz = Oxxz-Lambdaz;
FullOyyx = Oyyx-Lambdax;
FullOyyz = Oyyz-Lambdaz;
FullOzzx = Ozzx-Lambdax;
FullOzzy = Ozzy-Lambday;
FullOxxx = Oxxx-3.*Lambdax;
FullOyyy = Oyyy-3.*Lambday;
FullOzzz = Ozzz-3.*Lambdaz;

ScatO = ((k0.^8.*epsd.^2./(3780.*pi.*eps0.^2.*vd2.*mu0)).*(6.*abs(FullOxyz).^2+3.*abs(FullOxxy).^2+...  % Electric Octupole % электрический октуполь
    3.*abs(FullOxxz).^2+3.*abs(FullOyyx).^2+3.*abs(FullOyyz).^2+3.*abs(FullOzzx).^2+...
    3.*abs(FullOzzy).^2+abs(FullOxxx).^2+abs(FullOyyy).^2+abs(FullOzzz).^2)) ./ Iinc ;          
                
% Scattering cross-section as sum of multipole components
% Сечение рассеяние как сумма мультипольных компонентов
ScatCS = (ScatD + Scatm + ScatQ + ScatM + ScatO) ; 


lambda_nm = lambda .* norm_length;  % converting wavelenght from m to nm % переводим длину волны из метров в нанометры
H = H .* norm_length;               % converting height to nm % переводим высоту в нанометры


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPxTx = LocalMin( abs(TxK + Px), n_max );
maxPxTx = LocalMax( abs(TxK + Px), n_max );

fig1 = figure (1);

pl1 = @(n) plot (lambda_nm(n,:), abs(Px(n,:))./geomCS, ...
                 lambda_nm(n,:), abs(TxK(n,:) + Px(n,:))./geomCS, ...
                 lambda_nm(n,:), abs(TxK(n,:))./geomCS, ...
                 lambda_nm(n, maxPxTx{n}(:,1)), maxPxTx{n}(:,2)./geomCS, 'o', ...
                 lambda_nm(n, minPxTx{n}(:,1)), minPxTx{n}(:,2)./geomCS, 'o', ...
                 'LineWidth', LineWidth);
axis1 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([abs(Px)./geomCS, abs(TxK + Px)./geomCS, abs(TxK./geomCS)]))) ]);
xlab1 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab1 = @() ylabel ('dipole, a.u.','FontSize', FontSize);
leg1 =  @() legend( {'ED','ED+TD','TD'},'FontSize', FontSizeLeg);
tit1 =  @(n) title(strcat('Height h = ', num2str(H(n,1)), ' nm'),'FontSize', FontSize);
slider_toroidal( fig1, pl1, axis1, xlab1, ylab1, leg1, tit1, n, n_min, n_max, H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multipoles contributions to Extincntion
% Вклады мультиполей в экстинкцию
fig2 = figure(2);

pl2 = @(n) plot (lambda_nm(n,:), ExtPx(n,:)./geomCS, ...
                lambda_nm(n,:), ExtTx(n,:)./geomCS, ...
                lambda_nm(n,:), ExtQxz(n,:)./geomCS, ...
                lambda_nm(n,:), Extmy(n,:)./geomCS, ...
                lambda_nm(n,:), ExtMyz(n,:)./geomCS, ...
                lambda_nm(n,:), ExtOxzz(n,:)./geomCS, ...
                lambda_nm(n,:), ExtCS(n,:)./geomCS, ...
                'LineWidth', LineWidth);
axis2 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([ExtPx./geomCS, ExtTx./geomCS, ExtQxz./geomCS, Extmy./geomCS, ExtMyz./geomCS, ExtOxzz./geomCS, ExtCS./geomCS]))) ]);
tit2 = @(n) title(strcat('Multipoles contributions to Extincntion, h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab2 = @() xlabel ('Wavelength ,nm','FontSize', FontSize);
ylab2 = @() ylabel ('Multipoles contributions, um^2','FontSize', FontSize); % um - микрометры
leg2 = @() legend({'Px','Tx','Qxz','my','Myz','Oxzz','Total Ext-on Cross Sect. as Sum'},'FontSize', FontSizeLeg);
slider_toroidal( fig2, pl2, axis2, xlab2, ylab2, leg2, tit2, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All Cross-sections % Все сечения 
fig3 = figure(3);

pl3 = @(n) plot (lambda_nm(n,:), abs(ExtCS(n,:))./geomCS, ...
                lambda_nm(n,:), scat(n,:)./geomCS, ...
                lambda_nm(n,:), absCS(n,:)./geomCS, ...
                lambda_nm(n,:), (scat(n,:)+absCS(n,:))./geomCS, ...
                lambda_nm(n,:), ScatCS(n,:)./geomCS, ...
                'LineWidth', LineWidth);
axis3 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([abs(ExtCS)./geomCS, scat./geomCS, absCS./geomCS, (scat+absCS)./geomCS, ScatCS./geomCS]))) ]);
tit3 = @(n) title(strcat('Cross section, h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab3 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab3 = @() ylabel ('Cross section, mkm^2','FontSize', FontSize);
leg3 = @() legend({'Extincntion cross section','Scattering cross section from COMSOL','Absorption cross section from COMSOL', 'Extinction cross section from COMSOL (Abs+Scat)', 'Scatterning cross section'},'FontSize', FontSizeLeg );
slider_toroidal( fig3, pl3, axis3, xlab3, ylab3, leg3, tit3, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scattering cross-sections only % Только рассеяния
fig4 = figure (4);

pl4 = @(n) plot (lambda_nm(n,:), ScatCS(n,:)./geomCS, ...
                lambda_nm(n,:), scat(n,:)./geomCS, ...
                'LineWidth', LineWidth);
axis4 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([scat./geomCS, ScatCS./geomCS]))) ]);
tit4 = @(n) title(strcat('Cross section, h = ', num2str(H(n,1)),' nm' ),'FontSize', FontSize);
xlab4 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab4 = @() ylabel ('Cross section, um^2','FontSize', FontSize);
leg4 = @() legend({'Scattering cross section', 'Scattering cross section from COMSOL'},'FontSize', FontSizeLeg);
slider_toroidal( fig4, pl4, axis4, xlab4, ylab4, leg4, tit4, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig5 = figure (5);

pl5 = @(n) plot (lambda_nm(n,:), ScatD(n,:)./geomCS, ...
                lambda_nm(n,:), Scatm(n,:)./geomCS, ...
                lambda_nm(n,:), ScatQ(n,:)./geomCS, ...
                lambda_nm(n,:), ScatM(n,:)./geomCS, ...
                lambda_nm(n,:), ScatO(n,:)./geomCS, ...
                lambda_nm(n,:), ScatCS(n,:)./geomCS, ...
                'LineWidth', LineWidth);
axis5 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([ScatD./geomCS, Scatm./geomCS, ScatQ./geomCS, ScatM./geomCS, ScatO./geomCS, ScatCS./geomCS]))) ]);
tit5 = @(n) title(strcat('Multipoles Contributions to Scattering, h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab5 = @() xlabel ('Wavelenght, nm','FontSize', FontSize);
ylab5 = @() ylabel ('Multipoles Contributions, um^2','FontSize', FontSize);
leg5 = @() legend({'scat D ', 'scat m', 'scat Q', 'scat M', 'scat O', 'Sum Scat'},'FontSize', FontSizeLeg);
slider_toroidal( fig5, pl5, axis5, xlab5, ylab5, leg5, tit5, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% At the MaxValue 
% 1st line - maximum value                                        % 1 строка - максимальное значение величины
% 2nd line - frequency, corresponding to this value               % 2 строка - частота соответствующая этому значению
% 3rd line - wavelenght, corresponding to this value              % 3 строка - длина волны соответствующая этому значению
% 4th line - value of parameter, corresponding to maximum value   % 4 строка - значение параметра которому соответствует максимальное значение

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_lim = n_max;  % can be changed if you need to cut off all results after n_lim

Maxscat = MaxValue(scat./geomCS, fre, H, n_max); 

fig6 = figure(6);
set(fig6, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Scattering Cross-Section (Comsol) according to Height

subplot(1,2,1); hold on;
plot(Maxscat(4,n_min:n_lim), Maxscat(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max Scat. C-S (COMSOL) according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max Scat. Cross-section, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max Scattering Cross-Section (Comsol) according to Height
subplot(1,2,2); hold on;
plot(Maxscat(4,n_min:n_lim), Maxscat(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength corresponding to the Max Scat. C-S (COMSOL)','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Res Wavelength, nm','FontSize', FontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MaxAbs = MaxValue(absCS./geomCS, fre, H, n_max); 

fig7 = figure(7);
set(fig7, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Absorption Cross-Section (Comsol) according to Height

subplot(1,2,1); hold on;
plot(MaxAbs(4,n_min:n_lim), MaxAbs(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max Abs. C-S (COMSOL) according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max Abs. Cross-section, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max Absorption Cross-Section (Comsol) according to Height
subplot(1,2,2); hold on;
plot(MaxAbs(4,n_min:n_lim), MaxAbs(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength corresponding to the Max Abs. C-S (COMSOL)','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Res Wavelength, nm','FontSize', FontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MaxScatD = MaxValue(ScatD./geomCS, fre, H, n_max); 

fig8 = figure(8);
set(fig8, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Scattering component by Total Electric Dipole according to Height

subplot(1,2,1); hold on;
plot(MaxScatD(4,n_min:n_lim), MaxScatD(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max TED according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max TED, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max of scattering component by Total Electric Dipole according to Height
subplot(1,2,2); hold on;
plot(MaxScatD(4,n_min:n_lim), MaxScatD(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength corresponding to the Max TED','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Res Wavelength, nm','FontSize', FontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Maxm = MaxValue(Scatm./geomCS, fre, H, n_max); 

fig9 = figure(9);
set(fig9, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Scattering component by Magnetic Dipole according to Height
subplot(1,2,1); hold on;
plot(Maxm(4,n_min:n_lim), Maxm(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max magnetic dipole according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max magnetic dipole, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max Magnetic Dipole according to Height
subplot(1,2,2); hold on;
plot(Maxm(4,n_min:n_lim), Maxm(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength coresponding to the max magnetic dipole','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Wavelength, nm','FontSize', FontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxQ = MaxValue(ScatQ./geomCS, fre, H, n_max); 

fig10 = figure(10);
set(fig10, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Scattering component by Electric Quadrupole according to Height
subplot(1,2,1); hold on;
plot(MaxQ(4,n_min:n_lim), MaxQ(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max Electric Quadrupole according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max Electric Quadrupole, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max Electric Quadrupole according to Height
subplot(1,2,2); hold on;
plot(MaxQ(4,n_min:n_lim), MaxQ(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength coresponding to the Max Q','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Wavelength, nm','FontSize', FontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxM = MaxValue(ScatM./geomCS, fre, H, n_max); 

fig11 = figure(11);
set(fig11, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: max of Scattering component by Magnetic Quadrupole according to Height
subplot(1,2,1); hold on;
plot(MaxM(4,n_min:n_lim), MaxM(1,n_min:n_lim), 'r--o', ...
                'LineWidth', LineWidth); 
title('Max magnetic Quadrupole according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Max Magnetic Quadrupole, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to max Magnetic Quadrupole to Height
subplot(1,2,2); hold on;
plot(MaxM(4,n_min:n_lim), MaxM(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength coresponding to the Max M','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Wavelength, nm','FontSize', FontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MaxTxkPx = MaxValue(abs(TxK+Px), fre, H, n_max);
indexfre  = ones(1, n_max);
for i = 1 : n_max
    indexfre(1,i) = find(fre(1, :) == MaxTxkPx(2,i));
end

tempPx = ones(1, n_max);
for i = 1 : n_max
    tempPx(1,i) = Px(i, indexfre(1,i));
end

fig12 = figure(12);
set(fig12, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

% Dependence: (TD+ED)/ED in the point of maximum TD+ED according to Height
% Зависимость отношения суммы тороидального и электрического дипольного моментов в точке максимума TD+ED от высоты
subplot(1,2,1); hold on;
plot(MaxTxkPx(4,n_min:n_lim), MaxTxkPx(1,n_min:n_lim) ./ abs(tempPx(n_min:n_lim)), 'r--o', ...
                'LineWidth', LineWidth); 
title('(TD+ED)/ED at the Freq of Max (TD+ED) according to h','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('TD+ED/ED, a.u.','FontSize', FontSize);

% Dependence: Wavelenght, corresponding to (TD+ED)/ED in the point of maximum TD+ED to Height
% Зависимость длины волны, соответствующей отношению суммы тороидального и электрического дипольного моментов в точке максимума TD+ED, от высоты
subplot(1,2,2); hold on;
plot(MaxTxkPx(4,n_min:n_lim), MaxTxkPx(3,n_min:n_lim).*norm_length, 'r--o', ...
                'LineWidth', LineWidth);  
title('the Wavelength coresponding to the Max (TD+ED)/ED','FontSize', FontSize);
xlabel ('Height, nm','FontSize', FontSize);
ylabel ('Wavelength, nm','FontSize', FontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig15 = figure (15);

pl15 = @(n) plot (lambda_nm(n,:), ForScat(n,:).*1e16, ...                 
                 lambda_nm(n,:), ForScat(n,:)./BackScat(n,:), ...               
                 'LineWidth', LineWidth);
             %lambda_nm(n,:), BackScat(n,:).*1e16, ...    
             %lambda_nm(n,:), BackScat(n,:)./ForScat(n,:), ...
axis15 = @(n) axis([-inf, Inf, -inf, Inf]);
tit15 = @(n) title(strcat('Scattering (Patterns), h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab15 = @() xlabel ('Wavelenght, nm','FontSize', FontSize);
ylab15 = @() ylabel ('Far-field Scattering, a.u.', 'FontSize', FontSize);
leg15 = @() legend({'Forward Scattering', 'Forward/Backward Scattering', 'For/Back', 'Back/For'},'FontSize', FontSizeLeg);
slider_toroidal( fig15, pl15, axis15, xlab15, ylab15, leg15, tit15, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fig16 = figure (16); 

pl16 = @(n) plot (lambda_nm(n,:), ForScatPoint(n,:).*1e5, ...                 
                 lambda_nm(n,:), ForScatPoint(n,:)./BackScatPoint(n,:), ...               
                 'LineWidth', LineWidth);
             %lambda_nm(n,:), BackScat(n,:).*1e16, ...    
             %lambda_nm(n,:), BackScat(n,:)./ForScat(n,:), ...
axis16 = @(n) axis([-inf, Inf, -inf, Inf]);
tit16 = @(n) title(strcat('Scattering (Patterns) (Point), h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab16 = @() xlabel ('Wavelenght, nm','FontSize', FontSize);
ylab16 = @() ylabel ('Far-field Scattering, a.u.', 'FontSize', FontSize);
leg16 = @() legend({'Forward Scattering', 'Forward/Backward Scattering'},'FontSize', FontSizeLeg);
slider_toroidal( fig16, pl16, axis16, xlab16, ylab16, leg16, tit16, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig17 = figure (17); %phases

pl17 = @(n) plot (lambda_nm(n,:), angle(Dx(n,:))./pi, ...
                lambda_nm(n,:), angle(my(n,:))./pi, ...
                lambda_nm(n,:), angle(Dx(n,:))./pi+angle(my(n,:))./pi, ...
                'LineWidth', LineWidth);
            
%               lambda_nm(n,:), angle(Tx(n,:))./geomCS, ...
%               lambda_nm(n,:), angle(my(n,:))./geomCS, ...
%               lambda_nm(n,:), angle(my(n,:))./geomCS, ...
%               lambda_nm(n,:), angle(my(n,:))./geomCS, ...
axis17 = @(n) axis([-inf, Inf, -inf, inf]);
tit17 = @(n) title(strcat('Phases of moments, h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab17 = @() xlabel ('Wavelenght, nm','FontSize', FontSize);
ylab17 = @() ylabel ('Phase, pi','FontSize', FontSize);
leg17 = @() legend({'angle(Dx)', 'angle(my)','Sum'},'FontSize', FontSizeLeg);
slider_toroidal( fig17, pl17, axis17, xlab17, ylab17, leg17, tit17, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{

MaxTxk = MaxValue(abs(TxK)./geomCS, fre, H, n_max); % вывод нужных значений в файлы % writing values at files

scat_TED_to_H(:,1)=MaxScatD(4,n_min:n_lim);
scat_TED_to_H(:,2)=MaxScatD(1,n_min:n_lim);
scat_TED_to_H(:,3)=MaxScatD(3,n_min:n_lim);
dlmwrite('scat_TED_to_H_PP_100.csv', [scat_TED_to_H], 'delimiter', ',');

scat_m_to_H(:,1)=Maxm(4,n_min:n_lim);
scat_m_to_H(:,2)=Maxm(1,n_min:n_lim);
scat_m_to_H(:,3)=Maxm(3,n_min:n_lim);
dlmwrite('scat_m_to_H_PP_100.csv', [scat_m_to_H], 'delimiter', ',');

scat_TDEDRel_to_H(:,1)=MaxTxkPx(4,n_min:n_lim) ./ abs(tempPx(n_min:n_lim));
scat_TDEDRel_to_H(:,2)=MaxTxkPx(1,n_min:n_lim);
scat_TDEDRel_to_H(:,3)=MaxTxkPx(3,n_min:n_lim).*norm_length;
dlmwrite('scat_TDEDRel_to_H_PP_100.csv', [scat_TDEDRel_to_H], 'delimiter', ',');


scat_MQ_to_H(:,1)=MaxM(4,n_min:n_lim);
scat_MQ_to_H(:,2)=MaxM(1,n_min:n_lim);
scat_MQ_to_H(:,3)=MaxM(3,n_min:n_lim).*norm_length;
dlmwrite('scat_MQ_to_H_PP_100.csv', [scat_MQ_to_H], 'delimiter', ',');

scat_Q_to_H(:,1)=MaxQ(4,n_min:n_lim);
scat_Q_to_H(:,2)=MaxQ(1,n_min:n_lim);
scat_Q_to_H(:,3)=MaxQ(3,n_min:n_lim).*norm_length;
dlmwrite('scat_Q_to_H_PP_100.csv', [scat_Q_to_H], 'delimiter', ',');

scat_Tx_to_H(:,1)=MaxTxk(4,n_min:n_lim);
scat_Tx_to_H(:,2)=MaxTxk(1,n_min:n_lim);
scat_Tx_to_H(:,3)=MaxTxk(3,n_min:n_lim);
dlmwrite('scat_Tx_to_H_PP_100.csv', [scat_Tx_to_H], 'delimiter', ',');


%}