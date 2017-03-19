% Multipole Decomposotion with toroidal moment separation.
% ver 5.1
%19.03.2017
% added scattering patterns mechanism, to build scattering patterns via field distribution. Initially commented, use if needed.

clc
clear all;
%% Control panel for scattering patterns
% 0 to disable, any another number to activate
% Неплохо б сделать граф интерфейс в след версиях
Electric_dipole=1;
Magnetic_dipole=1;
Electric_quadrupole=1;
Magnetic_quadrupole=1;
Electrical_octupole=1;
Toroidal_dipole=1;
%%
n_min = 1;  % minimum graphic slider parameter and first considered parameter                   
            % минимальное значение слайдера; случаи с меньшим номером параметра не будут выводиться
f_min = 1;
n = 1;      % The original (first after fig creating) graphic slider parameter                  
            % первоначальное значение слайдера

FontSize = 16;          % font size at the titles % размер шрифта в подписях
FontSizeLeg = 9;        % font size at the legends % размер шрифта в легендах
LineWidth = 2.3;        % Line Width at the graphics % толщина линии на графиках            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Independent Variables
rad = 200e-9;   % base edge of pyramid/parallelepiped/etc or radius of sphere/cilinder/cone/etc, in METERS, MANUALLY, for Normalization! MUST BE CHECKED  
                % сторона основания пирамиды В МЕТРАХ, ПРОВЕРЯТЬ если нормируем на эффективное сечение!
%geomCS = rad/2; % effective geometrical cross-section for normalization, if it has square shape 
                % эффективное поперечное сечение рассеяния для нормировки на него
%geomCS = pi*rad^2; %effective geometrical cross-section for normalization, if it has circle shape. 
geomCS = 1e-18; % for usual normalization by nm  
%                 % нормировка на нм обычная!
% geomCS = 1e-12; % for usual normalization by um  
%                 % нормировка на мкм обычная!

%PS ЕСЛИ НАДО ДЕЛИТЬ НА КОНСТАНТУ РАССКОМЕНТИТЬ СТРОКИ В РАЙОНЕ 92-ой!
E0x = 1; 
E0 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
mu0 = 1.25663706 * 10^(-6); % Vacuum permeability (magnetic constant)     % магнитная постоянная 
eps0 = 8.85 * 10^(-12);     % Vacuum permittivity (electric constant)     % электрическая постоянная
c = 2.998e+8;               % speed of light                              % скорость света
% Z = 120*pi;               % for phase diagrams, still legacy parameter  % для фазовых диаграмм
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_length      = dlmread ('dim_value.dat');   % value to multiply for m -> nm converting  
                                                % величина, на которую домножаем, чтобы перевести метры в нанометры
dim_char         = fileread ('dim_name.dat');
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
Oxxy             = dlmread ('Oxxy.dat'); Oxyx=Oxxy; Oyxx=Oxxy;
Oxxz             = dlmread ('Oxxz.dat'); Oxzx=Oxxz; Ozxx=Oxxz;
Oyyx             = dlmread ('Oyyx.dat'); Oyxy=Oyyx; Oxyy=Oyyx;
Oyyy             = dlmread ('Oyyy.dat');
Oyyz             = dlmread ('Oyyz.dat'); Oyzy=Oyyz; Ozyy=Oyyz;
Ozzx             = dlmread ('Ozzx.dat'); Ozxz=Ozzx; Oxzz=Ozzx;
Ozzy             = dlmread ('Ozzy.dat'); Ozyz=Ozzy; Oyzz=Ozzy;
Ozzz             = dlmread ('Ozzz.dat');
Oxyz             = dlmread ('Oxyz.dat'); Oxzy=Oxyz; Oyxz=Oxyz; Oyzx=Oxyz; Ozxy=Oxyz; Ozyx=Oxyz;
Lambdax          = dlmread ('Lambdax.dat');
Lambday          = dlmread ('Lambday.dat');
Lambdaz          = dlmread ('Lambdaz.dat');
absCS            = dlmread ('absCS.dat');
scat             = dlmread ('scat.dat');
ForScat          = dlmread ('ForScat.dat');
BackScat         = dlmread ('BackScat.dat');
ForScatPoint     = dlmread ('ForScatPoint.dat');
BackScatPoint    = dlmread ('BackScatPoint.dat');

%geom_CS = geomCS.*H;  % effective geometrical cross-section for normalization, if it has square shape 
                      % эффективное поперечное сечение рассеяния для нормировки на него
geom_CS = geomCS.*(H.*0+1);  % РАССКОМЕНТИТЬ, ЕСЛИ СЕЧЕНИЕ НОРМИРОВКИ НЕ ЗАВИСИТ ОТ высоты и надо нормировать на geomCS!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_max = size(Px, 1);
f_max = size(Px, 2);
epsd = 1;   % Permittivity of environment outside of the particle 
            % диэлектричекая проницаемость среды снаружи частицы
lambda = c./fre;    % Wavelength as speed of light divided by frequency 
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
con=(((k0.^2).*exp(1i*kd))/(4*pi*eps0));% Change of notation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ud=1/sqrt(eps0*eps*mu0);                           % Light speed in homogeneous medium with eps
%% Extinction components for the wave, directed along Z and polarized along X
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

%% Extinction as sum of multipole components
% Экстинкция как сумма мультипольных составляющиx
ExtCS = ExtPx+ExtTx+ExtQxz+Extmy+ExtMyz+ExtOxzz;

%% Components of Total Electric Dipole Moment (TED)
% Компоненты полного электрического дипольного момента
Dx = Px - (1i.*k0.*epsd./c).*Tx;
Dy = Py - (1i.*k0.*epsd./c).*Ty;
Dz = Pz - (1i.*k0.*epsd./c).*Tz;

%% I incident - Poynting vector of incident wave, power flux density
% вектор пойнтинга падающей волны, плотность потока мощности
Iinc = ((eps0*epsd/mu0)^0.5*E0/2);

%% Contributions of separate miltipoles at the power scattering, normalized by Poynting vector
% Вклады отдельных мультиполей в рассеяние мощности нормированные на вектор пойнтинга
ScatD = (k0.^4./(12.*pi.*eps0.^2.*vd2.*mu0)) .* (abs(Dx).^2 + abs(Dy).^2 + abs(Dz).^2) ./ Iinc;     % Electric Dipole % диполь

Scatm = ((k0.^4.*epsd./(12.*pi*eps0.*vd2)) .* (abs(mx).^2+abs(my).^2+abs(mz).^2)) ./ Iinc;  % Magnetic Dipole % магнитный диполь

ScatQ=((k0.^6.*epsd./(1440.*pi.*eps0.^2.*vd2.*mu0))...                                  % Electric Quadrupole % электрический квадруполь
    .*(abs(Qxx).^2+abs(Qxy).^2+abs(Qxz).^2+abs(Qyx).^2+abs(Qyy).^2+abs(Qyz).^2+...
    abs(Qzx).^2+abs(Qzy).^2+abs(Qzz).^2)) ./ Iinc;                                          

ScatM = ((k0.^6.*epsd.^2./(160.*pi*eps0.^1.*vd2)) ...                                   % Magnetic quadrupole % магнитный квадруполь
    .*(abs(Mxx).^2+abs(Mxy).^2+abs(Mxz).^2+abs(Myx).^2+abs(Myy).^2+abs(Myz).^2+...
    abs(Mzx).^2+abs(Mzy).^2+abs(Mzz).^2) ) ./ Iinc;                                         

%% octupole parts from COMSOL % запчасти от октуполя из комсола
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
                
%% Scattering cross-section as sum of multipole components
% Сечение рассеяние как сумма мультипольных компонентов
ScatCS = (ScatD + Scatm + ScatQ + ScatM + ScatO) ; 


lambda_nm = lambda .* norm_length;  % converting Wavelength from m to nm % переводим длину волны из метров в нанометры
H = H .* norm_length;               % converting height to nm % переводим высоту в нанометры

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
minPxTx = LocalMin( abs(TxK + Px), n_max );
maxPxTx = LocalMax( abs(TxK + Px), n_max );


fig1 = figure (1);

pl1 = @(n) plot (lambda_nm(n,:), abs(Px(n,:))./geom_CS(n,:), ...
                 lambda_nm(n,:), abs(TxK(n,:) + Px(n,:))./geom_CS(n,:), ...
                 lambda_nm(n,:), abs(TxK(n,:))./geom_CS(n,:), ...
                 lambda_nm(n, maxPxTx{n}(:,1)), maxPxTx{n}(:,2)./geom_CS(n,1), 'o', ...
                 lambda_nm(n, minPxTx{n}(:,1)), minPxTx{n}(:,2)./geom_CS(n,1), 'o', ...
                 'LineWidth', LineWidth);
axis1 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([abs(Px)./min(geom_CS(:,1)), abs(TxK + Px)./min(geom_CS(:,1)), abs(TxK./min(geom_CS(:,1)))]))) ]);
xlab1 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab1 = @() ylabel ('dipole, a.u.','FontSize', FontSize);
leg1 =  @() legend( {'ED','ED+TD','TD'},'FontSize', FontSizeLeg);
tit1 =  @(n) title (strcat('Height h = ',32, num2str(H(n,1)), 32, dim_char),'FontSize', FontSize);
slider_toroidal( fig1, pl1, axis1, xlab1, ylab1, leg1, tit1, n, n_min, n_max, H);
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Multipoles contributions to Extincntion
% Вклады мультиполей в экстинкцию
fig2 = figure(2);

pl2 = @(n) plot (lambda_nm(n,:), ExtPx(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ExtTx(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ExtQxz(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), Extmy(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ExtMyz(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ExtOxzz(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ExtCS(n,:)./geom_CS(n,:), ...
                'LineWidth', LineWidth);
axis2 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([ExtPx./min(geom_CS(:,1)), ExtTx./min(geom_CS(:,1)), ExtQxz./min(geom_CS(:,1)), Extmy./min(geom_CS(:,1)), ExtMyz./min(geom_CS(:,1)), ExtOxzz./min(geom_CS(:,1)), ExtCS./min(geom_CS(:,1))]))) ]);
tit2 = @(n) title (strcat('Multipoles contributions to Extincntion, h = ',32, num2str(H(n,1)), 32, dim_char),'FontSize', FontSize);
xlab2 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab2 = @() ylabel ('Multipoles contributions, nm^2','FontSize', FontSize); % um - микрометры
leg2 = @() legend({'Px','Tx','Qxz','my','Myz','Oxzz','Total Ext-on Cross Sect. as Sum'},'FontSize', FontSizeLeg);
slider_toroidal( fig2, pl2, axis2, xlab2, ylab2, leg2, tit2, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All Cross-sections % Все сечения 
fig3 = figure(3);

pl3 = @(n) plot (lambda_nm(n,:), abs(ExtCS(n,:))./geom_CS(n,:), ...
                lambda_nm(n,:), scat(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), absCS(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), (scat(n,:)+absCS(n,:))./geom_CS(n,:), ...
                lambda_nm(n,:), ScatCS(n,:)./geom_CS(n,:), ...
                'LineWidth', LineWidth);
axis3 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([abs(ExtCS)./min(geom_CS(:,1)), scat./min(geom_CS(:,1)), absCS./min(geom_CS(:,1)), (scat+absCS)./min(geom_CS(:,1)), ScatCS./min(geom_CS(:,1))]))) ]);
tit3 = @(n) title (strcat('Cross section, h = ',32, num2str(H(n,1)), 32, dim_char),'FontSize', FontSize);
xlab3 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab3 = @() ylabel ('Cross section, nm^2','FontSize', FontSize);
leg3 = @() legend({'Extincntion cross section','Scattering cross section from COMSOL','Absorption cross section from COMSOL', 'Extinction cross section from COMSOL (Abs+Scat)', 'Scatterning cross section'},'FontSize', FontSizeLeg );
slider_toroidal( fig3, pl3, axis3, xlab3, ylab3, leg3, tit3, n, n_min, n_max, H );
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scattering cross-sections only % Только рассеяния
fig4 = figure (4);

pl4 = @(n) plot (lambda_nm(n,:), ScatCS(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), scat(n,:)./geom_CS(n,:), ...
                'LineWidth', LineWidth);
axis4 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([scat./min(geom_CS(:,1)), ScatCS./min(geom_CS(:,1))]))) ]);
tit4 = @(n) title (strcat('Cross section, h = ',32, num2str(H(n,1)), 32, dim_char),'FontSize', FontSize);
xlab4 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab4 = @() ylabel ('Cross section, nm^2','FontSize', FontSize);
leg4 = @() legend({'Scattering cross section', 'Scattering cross section from COMSOL'},'FontSize', FontSizeLeg);
slider_toroidal( fig4, pl4, axis4, xlab4, ylab4, leg4, tit4, n, n_min, n_max, H );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig5 = figure (5);

pl5 = @(n) plot (lambda_nm(n,:), ScatD(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), Scatm(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ScatQ(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ScatM(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ScatO(n,:)./geom_CS(n,:), ...
                lambda_nm(n,:), ScatCS(n,:)./geom_CS(n,:), ...
                'LineWidth', LineWidth);
axis5 = @(n) axis([-inf, Inf, -inf, 1.05*(max(max([ScatD./min(geom_CS(:,1)), Scatm./min(geom_CS(:,1)), ScatQ./min(geom_CS(:,1)), ScatM./min(geom_CS(:,1)), ScatO./min(geom_CS(:,1)), ScatCS./min(geom_CS(:,1))]))) ]);
tit5 = @(n) title (strcat('Multipoles Contributions to Scattering, h = ',32, num2str(H(n,1)), 32, dim_char),'FontSize', FontSize);
xlab5 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab5 = @() ylabel ('Multipoles Contributions, nm^2','FontSize', FontSize);
leg5 = @() legend({'scat D ', 'scat m', 'scat Q', 'scat M', 'scat O', 'Sum Scat'},'FontSize', FontSizeLeg);
slider_toroidal( fig5, pl5, axis5, xlab5, ylab5, leg5, tit5, n, n_min, n_max, H );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % At the MaxValue 
% % 1st line - maximum value                                        % 1 строка - максимальное значение величины
% % 2nd line - frequency, corresponding to this value               % 2 строка - частота соответствующая этому значению
% % 3rd line - Wavelength, corresponding to this value              % 3 строка - длина волны соответствующая этому значению
% % 4th line - value of parameter, corresponding to maximum value   % 4 строка - значение параметра которому соответствует максимальное значение
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_lim = n_max;  % can be changed if you need to cut off all results after n_lim
% 
% Maxscat = MaxValue(scat./geom_CS(:,:), fre, H, n_max); 
% 
% fig6 = figure(6);
% set(fig6, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Scattering Cross-Section (Comsol) according to Height
% 
% subplot(1,2,1); hold on;
% plot(Maxscat(4,n_min:n_lim), Maxscat(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title ('Max Scat. C-S (COMSOL) according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max Scat. Cross-section, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max Scattering Cross-Section (Comsol) according to Height
% subplot(1,2,2); hold on;
% plot(Maxscat(4,n_min:n_lim), Maxscat(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength corresponding to the Max Scat. C-S (COMSOL)','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Res Wavelength, nm','FontSize', FontSize);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% MaxAbs = MaxValue(absCS./geom_CS(:,:), fre, H, n_max); 
% 
% fig7 = figure(7);
% set(fig7, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Absorption Cross-Section (Comsol) according to Height
% 
% subplot(1,2,1); hold on;
% plot(MaxAbs(4,n_min:n_lim), MaxAbs(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('Max Abs. C-S (COMSOL) according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max Abs. Cross-section, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max Absorption Cross-Section (Comsol) according to Height
% subplot(1,2,2); hold on;
% plot(MaxAbs(4,n_min:n_lim), MaxAbs(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength corresponding to the Max Abs. C-S (COMSOL)','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Res Wavelength, nm','FontSize', FontSize);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MaxScatD = MaxValue(ScatD./geom_CS(:,:), fre, H, n_max); 
% 
% fig8 = figure(8);
% set(fig8, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Scattering component by Total Electric Dipole according to Height
% 
% subplot(1,2,1); hold on;
% plot(MaxScatD(4,n_min:n_lim), MaxScatD(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('Max TED according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max TED, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max of scattering component by Total Electric Dipole according to Height
% subplot(1,2,2); hold on;
% plot(MaxScatD(4,n_min:n_lim), MaxScatD(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength corresponding to the Max TED','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Res Wavelength, nm','FontSize', FontSize);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maxm = MaxValue(Scatm./geom_CS(:,:), fre, H, n_max); 
% 
% fig9 = figure(9);
% set(fig9, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Scattering component by Magnetic Dipole according to Height
% subplot(1,2,1); hold on;
% plot(Maxm(4,n_min:n_lim), Maxm(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('Max magnetic dipole according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max magnetic dipole, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max Magnetic Dipole according to Height
% subplot(1,2,2); hold on;
% plot(Maxm(4,n_min:n_lim), Maxm(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength coresponding to the max magnetic dipole','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Wavelength, nm','FontSize', FontSize);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaxQ = MaxValue(ScatQ./geom_CS(:,:), fre, H, n_max); 
% 
% fig10 = figure(10);
% set(fig10, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Scattering component by Electric Quadrupole according to Height
% subplot(1,2,1); hold on;
% plot(MaxQ(4,n_min:n_lim), MaxQ(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('Max Electric Quadrupole according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max Electric Quadrupole, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max Electric Quadrupole according to Height
% subplot(1,2,2); hold on;
% plot(MaxQ(4,n_min:n_lim), MaxQ(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength coresponding to the Max Q','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Wavelength, nm','FontSize', FontSize);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaxM = MaxValue(ScatM./geom_CS(:,:), fre, H, n_max); 
% 
% fig11 = figure(11);
% set(fig11, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: max of Scattering component by Magnetic Quadrupole according to Height
% subplot(1,2,1); hold on;
% plot(MaxM(4,n_min:n_lim), MaxM(1,n_min:n_lim), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('Max magnetic Quadrupole according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Max Magnetic Quadrupole, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to max Magnetic Quadrupole to Height
% subplot(1,2,2); hold on;
% plot(MaxM(4,n_min:n_lim), MaxM(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength coresponding to the Max M','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Wavelength, nm','FontSize', FontSize);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% MaxTxkPx = MaxValue(abs(TxK+Px), fre, H, n_max);
% indexfre  = ones(1, n_max);
% for i = 1 : n_max
%     indexfre(1,i) = find(fre(1, :) == MaxTxkPx(2,i));
% end
% 
% tempPx = ones(1, n_max);
% for i = 1 : n_max
%     tempPx(1,i) = Px(i, indexfre(1,i));
% end
% 
% fig12 = figure(12);
% set(fig12, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);
% 
% % Dependence: (TD+ED)/ED in the point of maximum TD+ED according to Height
% % Зависимость отношения суммы тороидального и электрического дипольного моментов в точке максимума TD+ED от высоты
% subplot(1,2,1); hold on;
% plot(MaxTxkPx(4,n_min:n_lim), MaxTxkPx(1,n_min:n_lim) ./ abs(tempPx(n_min:n_lim)), 'r--o', ...
%                 'LineWidth', LineWidth); 
% title('(TD+ED)/ED at the Freq of Max (TD+ED) according to h','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('TD+ED/ED, a.u.','FontSize', FontSize);
% 
% % Dependence: Wavelength, corresponding to (TD+ED)/ED in the point of maximum TD+ED to Height
% % Зависимость длины волны, соответствующей отношению суммы тороидального и электрического дипольного моментов в точке максимума TD+ED, от высоты
% subplot(1,2,2); hold on;
% plot(MaxTxkPx(4,n_min:n_lim), MaxTxkPx(3,n_min:n_lim).*norm_length, 'r--o', ...
%                 'LineWidth', LineWidth);  
% title('the Wavelength coresponding to the Max (TD+ED)/ED','FontSize', FontSize);
% xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylabel ('Wavelength, nm','FontSize', FontSize);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fig15 = figure (15);
% 
% pl15 = @(n) plot (lambda_nm(n,:), ForScat(n,:)./geom_CS(n,:).*1e-5, ...                 
%                  lambda_nm(n,:), ForScat(n,:)./BackScat(n,:), ...               
%                  'LineWidth', LineWidth);
%              %lambda_nm(n,:), BackScat(n,:).*1e16, ...    
%              %lambda_nm(n,:), BackScat(n,:)./ForScat(n,:), ...
% axis15 = @(n) axis([-inf, Inf, -inf, Inf]);
% tit15 = @(n) title(strcat('Scattering (Patterns), h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
% xlab15 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
% ylab15 = @() ylabel ('Far-field Scattering, a.u.', 'FontSize', FontSize);
% leg15 = @() legend({'Forward Scattering', 'Forward/Backward Scattering'},'FontSize', FontSizeLeg);
% slider_toroidal( fig15, pl15, axis15, xlab15, ylab15, leg15, tit15, n, n_min, n_max, H );
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%{
fig17 = figure (17); %phases

pl17 = @(n) plot (lambda_nm(n,:), angle(-Dx(n,:))./pi, ...
                lambda_nm(n,:), angle(my(n,:))./pi, ...
                lambda_nm(n,:), angle(Qxz(n,:))./pi, ...
                'LineWidth', LineWidth);
            
%               lambda_nm(n,:), angle(Tx(n,:))./geom_CS(n,:), ...
%               lambda_nm(n,:), angle(my(n,:))./geom_CS(n,:), ...
%               lambda_nm(n,:), angle(my(n,:))./geom_CS(n,:), ...
%               lambda_nm(n,:), angle(my(n,:))./geom_CS(n,:), ...
axis17 = @(n) axis([-inf, Inf, -inf, inf]);
tit17 = @(n) title(strcat('Phases of moments, h = ', num2str(H(n,1)), ' nm' ),'FontSize', FontSize);
xlab17 = @() xlabel (strcat('Wavelength,', 32, dim_char),'FontSize', FontSize);
ylab17 = @() ylabel ('Phase, pi','FontSize', FontSize);
leg17 = @() legend({'TED', 'MD', 'Qxz', 'Qzx'},'FontSize', FontSizeLeg);
slider_toroidal( fig17, pl17, axis17, xlab17, ylab17, leg17, tit17, n, n_min, n_max, H );
%}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %{
% 
% MaxTxk = MaxValue(abs(TxK)./geom_CS(n,:), fre, H, n_max); % вывод нужных значений в файлы % writing values at files
% 
% scat_TED_to_H(:,1)=MaxScatD(4,n_min:n_lim);
% scat_TED_to_H(:,2)=MaxScatD(1,n_min:n_lim);
% scat_TED_to_H(:,3)=MaxScatD(3,n_min:n_lim);
% dlmwrite('scat_TED_to_H_PP_100.csv', [scat_TED_to_H], 'delimiter', ',');
% 
% scat_m_to_H(:,1)=Maxm(4,n_min:n_lim);
% scat_m_to_H(:,2)=Maxm(1,n_min:n_lim);
% scat_m_to_H(:,3)=Maxm(3,n_min:n_lim);
% dlmwrite('scat_m_to_H_PP_100.csv', [scat_m_to_H], 'delimiter', ',');
% 
% scat_TDEDRel_to_H(:,1)=MaxTxkPx(4,n_min:n_lim) ./ abs(tempPx(n_min:n_lim));
% scat_TDEDRel_to_H(:,2)=MaxTxkPx(1,n_min:n_lim);
% scat_TDEDRel_to_H(:,3)=MaxTxkPx(3,n_min:n_lim).*norm_length;
% dlmwrite('scat_TDEDRel_to_H_PP_100.csv', [scat_TDEDRel_to_H], 'delimiter', ',');
% 
% 
% scat_MQ_to_H(:,1)=MaxM(4,n_min:n_lim);
% scat_MQ_to_H(:,2)=MaxM(1,n_min:n_lim);
% scat_MQ_to_H(:,3)=MaxM(3,n_min:n_lim).*norm_length;
% dlmwrite('scat_MQ_to_H_PP_100.csv', [scat_MQ_to_H], 'delimiter', ',');
% 
% scat_Q_to_H(:,1)=MaxQ(4,n_min:n_lim);
% scat_Q_to_H(:,2)=MaxQ(1,n_min:n_lim);
% scat_Q_to_H(:,3)=MaxQ(3,n_min:n_lim).*norm_length;
% dlmwrite('scat_Q_to_H_PP_100.csv', [scat_Q_to_H], 'delimiter', ',');
% 
% scat_Tx_to_H(:,1)=MaxTxk(4,n_min:n_lim);
% scat_Tx_to_H(:,2)=MaxTxk(1,n_min:n_lim);
% scat_Tx_to_H(:,3)=MaxTxk(3,n_min:n_lim);
% dlmwrite('scat_Tx_to_H_PP_100.csv', [scat_Tx_to_H], 'delimiter', ',');
% 
% 
% 
%{
%% Scattering of power
% Диаграммы рассеяния мощности
if Electric_dipole==0
    Px=zeros(size(Px)); Py=Px; Pz=Px;  
end
if Magnetic_dipole==0
    mx=zeros(size(mx)); my=mx; mz=mx; 
end
if Electric_quadrupole==0
    Qxx=zeros(size(Qxx)); Qyy=Qxx; Qzz=Qxx; Qxy=Qxx; Qxz=Qxx; Qyx=Qxx; Qyz=Qxx; Qzx=Qxx; Qzy=Qxx;
end
if Magnetic_quadrupole==0
    Mxx=zeros(size(Mxx)); Myy=Mxx; Mzz=Mxx; Mxy=Mxx; Mxz=Mxx; Myx=Mxx; Myz=Mxx; Mzx=Mxx; Mzy=Mxx;
end
if Electrical_octupole==0
    Oxxx=zeros(size(Oxxx)); Oxxy=Oxxx; Oxxz=Oxxx; Oyyx=Oxxx; Oyyy=Oxxx; Oyyz=Oxxx; 
    Ozzx=Oxxx; Ozzy=Oxxx; Ozzz=Oxxx; Oxyz=Oxxx; 
    Oxyx=Oxxy; Oyxx=Oxxy; Oxzx=Oxxz; Ozxx=Oxxz; Oyxy=Oyyx; Oxyy=Oyyx; Oyzy=Oyyz; Ozyy=Oyyz; 
    Ozxz=Ozzx; Oxzz=Ozzx; Ozyz=Ozzy; Oyzz=Ozzy; Oxzy=Oxyz; Oyxz=Oxyz; Oyzx=Oxyz; Ozxy=Oxyz; Ozyx=Oxyz;
end
if Toroidal_dipole==0
    Tx=zeros(size(Tx)); Ty=Tx; Tz=Tx; 
end
%% Zeros (Attention! Depends on step of loops)
Ep_x=zeros(101,101); Ep_y=Ep_x; Ep_z=Ep_x; Mp_x=Ep_x; Mp_y=Ep_x; Mp_z=Ep_x;
EQ_x=Ep_x; EQ_y=Ep_x; EQ_z=Ep_x; MQ_x=Ep_x; MQ_y=Ep_x; MQ_z=Ep_x; EO_x=Ep_x; EO_y=Ep_x; EO_z=Ep_x; 
T_x=Ep_x; T_y=Ep_x; T_z=Ep_x; Field_x=Ep_x; Field_y=Ep_x; Field_z=Ep_x; P=Ep_x; x=Ep_x; y=Ep_x; z=Ep_x;
%% Calculations
del=pi/50;                      % Step of loops
for l=1:n_max
for f=1:f_max
j=0;                            % Counter of angle phi (rows), placed before 'phi' loop
for phi=0:del:2*pi              % Loop on rows
    j=j+1;                      % Next row
    k=0;                        % Counter of angle theta (columns), placed before 'theta' loop 
for    theta=0:del:2*pi         % Loop on columns
    k=k+1;                      % Next column
    nx=sin(theta)*cos(phi);     % x component of unit vector in spherical coordinates
    ny=sin(theta)*sin(phi);     % y component of unit vector in spherical coordinates
    nz=cos(theta);              % z component of unit vector in spherical coordinates
%% Electric Dipole    
    vec_Ep1=vector_product(Px(l,f),Py(l,f),Pz(l,f),nx,ny,nz);                          
    vec_Ep2=vector_product(nx,ny,nz,vec_Ep1(1),vec_Ep1(2),vec_Ep1(3));
Ep_x(j,k,l,f)=con(l,f)*vec_Ep2(1);
Ep_y(j,k,l,f)=con(l,f)*vec_Ep2(2);
Ep_z(j,k,l,f)=con(l,f)*vec_Ep2(3);
%% Magnetic Dipole
    vec_Mp1=vector_product(mx(l,f),my(l,f),mz(l,f),nx,ny,nz);
Mp_x(j,k,l,f)=(con(l,f)/ud)*vec_Mp1(1);
Mp_y(j,k,l,f)=(con(l,f)/ud)*vec_Mp1(2);
Mp_z(j,k,l,f)=(con(l,f)/ud)*vec_Mp1(3);
%% Electric Quadrupole
vec_EQ1(1)=Qxx(l,f)*nx+Qxy(l,f)*ny+Qxz(l,f)*nz;
vec_EQ1(2)=Qyx(l,f)*nx+Qyy(l,f)*ny+Qyz(l,f)*nz;
vec_EQ1(3)=Qzx(l,f)*nx+Qzy(l,f)*ny+Qzz(l,f)*nz;
    vec_EQ2= vector_product(nx,ny,nz,vec_EQ1(1),vec_EQ1(2),vec_EQ1(3));
    vec_EQ3= vector_product(nx,ny,nz,vec_EQ2(1),vec_EQ2(2),vec_EQ2(3));
EQ_x(j,k,l,f)=con(l,f)*(1i*kd(l,f)/6)*vec_EQ3(1);
EQ_y(j,k,l,f)=con(l,f)*(1i*kd(l,f)/6)*vec_EQ3(2);
EQ_z(j,k,l,f)=con(l,f)*(1i*kd(l,f)/6)*vec_EQ3(3);
%% Magnetic Quadrupole
vec_MQ1(1)=Mxx(l,f)*nx+Mxy(l,f)*ny+Mxz(l,f)*nz;
vec_MQ1(2)=Myx(l,f)*nx+Myy(l,f)*ny+Myz(l,f)*nz;
vec_MQ1(3)=Mzx(l,f)*nx+Mzy(l,f)*ny+Mzz(l,f)*nz;
    vec_MQ2= vector_product(nx,ny,nz,vec_MQ1(1),vec_MQ1(2),vec_MQ1(3));
MQ_x(j,k,l,f)=con(l,f)*((1i*kd(l,f))/(2*ud))*vec_MQ2(1);
MQ_y(j,k,l,f)=con(l,f)*((1i*kd(l,f))/(2*ud))*vec_MQ2(2);
MQ_z(j,k,l,f)=con(l,f)*((1i*kd(l,f))/(2*ud))*vec_MQ2(3);
%% Electric Octupole
O(1,1)=Oxxx(l,f)*nx+Oxxy(l,f)*ny+Oxxz(l,f)*nz;
O(1,2)=Oxyx(l,f)*nx+Oxyy(l,f)*ny+Oxyz(l,f)*nz;
O(1,3)=Oxzx(l,f)*nx+Oxzy(l,f)*ny+Oxzz(l,f)*nz;
O(2,1)=Oyxx(l,f)*nx+Oyxy(l,f)*ny+Oyxz(l,f)*nz;
O(2,2)=Oyyx(l,f)*nx+Oyyy(l,f)*ny+Oyyz(l,f)*nz;
O(2,3)=Oyzx(l,f)*nx+Oyzy(l,f)*ny+Oyzz(l,f)*nz;
O(3,1)=Ozxx(l,f)*nx+Ozxy(l,f)*ny+Ozxz(l,f)*nz;
O(3,2)=Ozyx(l,f)*nx+Ozyy(l,f)*ny+Ozyz(l,f)*nz;
O(3,3)=Ozzx(l,f)*nx+Ozzy(l,f)*ny+Ozzz(l,f)*nz;
    vec_EO2=dot_product_matrix_vector(O,nx,ny,nz);
    vec_EO3= vector_product(nx,ny,nz,vec_EO2(1),vec_EO2(2),vec_EO2(3));
    vec_EO4= vector_product(nx,ny,nz,vec_EO3(1),vec_EO3(2),vec_EO3(3));
EO_x(j,k,l,f)=con(l,f)*((kd(l,f)^2)/6)*vec_EO4(1);
EO_y(j,k,l,f)=con(l,f)*((kd(l,f)^2)/6)*vec_EO4(2);
EO_z(j,k,l,f)=con(l,f)*((kd(l,f)^2)/6)*vec_EO4(3);
%% Toroid moment
    vec_T1=vector_product(Tx(l,f),Ty(l,f),Tz(l,f),nx,ny,nz);                          
    vec_T2=vector_product(nx,ny,nz,vec_T1(1),vec_T1(2),vec_T1(3));
T_x(j,k,l,f)=con(l,f)*(1i*kd(l,f)/ud)*vec_T2(1);
T_y(j,k,l,f)=con(l,f)*(1i*kd(l,f)/ud)*vec_T2(2);
T_z(j,k,l,f)=con(l,f)*(1i*kd(l,f)/ud)*vec_T2(3);
%% Scattering Power
Field_x(j,k,l,f)=Ep_x(j,k,l,f)+Mp_x(j,k,l,f)+EQ_x(j,k,l,f)+MQ_x(j,k,l,f)+EO_x(j,k,l,f)+T_x(j,k,l,f);
Field_y(j,k,l,f)=Ep_y(j,k,l,f)+Mp_y(j,k,l,f)+EQ_y(j,k,l,f)+MQ_y(j,k,l,f)+EO_y(j,k,l,f)+T_y(j,k,l,f);
Field_z(j,k,l,f)=Ep_z(j,k,l,f)+Mp_z(j,k,l,f)+EQ_z(j,k,l,f)+MQ_z(j,k,l,f)+EO_z(j,k,l,f)+T_z(j,k,l,f);
% Field_x(j,k,l,f)=Ep_x(j,k,l,f);
% Field_y(j,k,l,f)=Ep_y(j,k,l,f);
% Field_z(j,k,l,f)=Ep_z(j,k,l,f);
P(j,k,l,f)=(1/2)*sqrt(eps0/mu0)*(((abs(Field_x(j,k,l,f))).^2)+((abs(Field_y(j,k,l,f))).^2)+(((abs(Field_z(j,k,l,f))).^2)));
%% Spherical Coordinates
x(j,k,l,f)=P(j,k,l,f).*(sin(theta).*cos(phi));
y(j,k,l,f)=P(j,k,l,f).*(sin(theta).*sin(phi));
z(j,k,l,f)=P(j,k,l,f).*cos(theta);
end
end
end
end

bound=@(n,f)[max(x(:,:,n,f)) max(y(:,:,n,f)) max(z(:,:,n,f)) abs(min(x(:,:,n,f))) abs(min(y(:,:,n,f))) abs(min(z(:,:,n,f)))];

figScat=figure;
phi=0:del:2*pi;
sir=@(n,f) surf(x(:,:,n,f),y(:,:,n,f),z(:,:,n,f)); % 3D график
pol1=@(n,f) polar(phi,P(1,:,n,f));                 % 'Cross-section in x-z surface' в полярных
pol2=@(n,f) polar(phi,P(26,:,n,f));                % 'Cross-section in y-z surface' в полярных
ss=mytitle(Electric_dipole, Magnetic_dipole,Electric_quadrupole,Magnetic_quadrupole,Electrical_octupole,Toroidal_dipole);
slider_scattering(figScat, sir, pol1, pol2, bound, ss, n, n_min, n_max, f, f_min, f_max, fre, H, norm_length)

%}
