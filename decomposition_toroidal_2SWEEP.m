% Мультипольная декомпозиция с выделением тороидального момента
% ver 1.1.3
% GitHub is wizard's hat

clc
clear all;

norm_length = 1e9; % величина, на которую домножаем, чтобы перевести метры в нанометры
% norm_length = 1e6; % величина, на которую домножаем, чтобы перевести метры в микрометры

n_max = 2;	% задавать вручную, число шагов в parametric sweep, максимальное значение слайдера
n_min = 1;	% минимальное значение слайдера
n = 1;     	% первоначальное значение слайдера


mu0 = 1.25663706 * 10^(-6);	% магнитная постоянная 
eps0 = 8.85 * 10^(-12);		% электрическая постоянная
c = 2.998e+8; 				% скорость света
% Z = 120*pi; 				% для фазовых диаграм

epsilon_tbl = dlmread ('Px.txt', '' ,5,0); 	% считывет данные из файла, обрезает 5 верхних строк
fre = epsilon_tbl(n:n_max:end,3); 		   	% считывает из 3 колонки начиная с n-ного элемента каждый n_max элемент
fre = repmat(fre,1,n_max);  % делаем из вектор-стобца fre матрицу, состоящую из n_max вектор-стобцов fre
length_fre = size(fre,1); 	% число шагов по частоте, длина вектора частот
lambda = c./fre;  			% длина волны - скорость света делить на частоту
epsd = 1;		 	% диэлектричекая проницаемость среды снаружи частицы
k0 = - 2*pi*fre/c; 	% волновой вектор 
				  	% минус потому что в COMSOL задан минус в мнимой экспоненте e^(-i k*x), описывающей падающее излучение

kd = k0 .* sqrt(epsd); % волновой вектор в среде
vd2 = c ./ sqrt(epsd); % скорость света в среде, vd для сечения рассеяния без зависимости от минуса в k0
% vd = c ./ sqrt(epsd);
vd = 2.*pi.*fre./( k0 .* sqrt(epsd) ); 	% тоже скорость света в среде. С этим vd хорошо считается экстинкция.

E0x = 1 ; 
E0 = 1;
rad = 100e-9;	% сторона основания пирамиды В МЕТРАХ, ПРОВЕРЯТЬ если нормируем на эффективное сечение!
geomCS = rad^2; % эффективное поперечное сечение рассеяния для нормировки на него
% geomCS = 1e-18; % нормировка на нм обычная!

% Height в метрах
H = ones(length_fre, n_max); 
for i = 1:1:n_max
H(:,i) = H(:,i) .* epsilon_tbl(i:n_max:end,2);
end

Px = ones(length_fre, n_max);
for i = 1:1:n_max
Px(:,i) = Px(:,i) .* epsilon_tbl(i:n_max:end,4);	% поляризация - плотность дипольного момента, X компонетна
end

epsilon_tbl = dlmread ('Py.txt', '' ,5,0);
Py = ones(length_fre, n_max);
for i = 1:1:n_max
Py(:,i) = Py(:,i) .* epsilon_tbl(i:n_max:end,4);	% считывает из 4 колонки начиная с n-ного элемента каждый n_max элемент
end

epsilon_tbl = dlmread ('Pz.txt', '' ,5,0);
Pz  = ones(length_fre, n_max);
for i = 1:1:n_max
Pz(:,i) = Pz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Tx.txt', '' ,5,0);
Tx  = ones(length_fre, n_max);
for i = 1:1:n_max
Tx(:,i) = Tx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Ty.txt', '' ,5,0);
Ty  = ones(length_fre, n_max);
for i = 1:1:n_max
Ty(:,i) = Ty(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Tz.txt', '' ,5,0);
Tz  = ones(length_fre, n_max);
for i = 1:1:n_max
Tz(:,i) = Tz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

% нормировка сечения рассеяния на вектор пойнтинга падающего излучения
epsilon_tbl=dlmread ('scat.txt', '' ,5,0);
scat  = ones(length_fre, n_max);
for i = 1:1:n_max
scat(:,i) = scat(:,i) .* epsilon_tbl(i:n_max:end,4) ./ 1.3E-3;	% сечение рассеяния в Ваттах
end                   % нормировка сечения рассеяния на вектор пойнтинга падающего излучения

epsilon_tbl=dlmread ('mx.txt', '' ,5,0);
mx  = ones(length_fre, n_max);
for i = 1:1:n_max
mx(:,i) = mx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('my.txt', '' ,5,0);
my  = ones(length_fre, n_max);
for i = 1:1:n_max
my(:,i) = my(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('mz.txt', '' ,5,0);
mz  = ones(length_fre, n_max);
for i = 1:1:n_max
mz(:,i) = mz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Mxx.txt', '' ,5,0);
Mxx  = ones(length_fre, n_max);
for i = 1:1:n_max
Mxx(:,i) = Mxx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Mxy.txt', '' ,5,0);
Mxy  = ones(length_fre, n_max);
for i = 1:1:n_max
Mxy(:,i) = Mxy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Mxz.txt', '' ,5,0);
Mxz  = ones(length_fre, n_max);
for i = 1:1:n_max
Mxz(:,i) = Mxz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Myx.txt', '' ,5,0);
Myx  = ones(length_fre, n_max);
for i = 1:1:n_max
Myx(:,i) = Myx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Myy.txt', '' ,5,0);
Myy  = ones(length_fre, n_max);
for i = 1:1:n_max
Myy(:,i) = Myy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Myz.txt', '' ,5,0);
Myz  = ones(length_fre, n_max);
for i = 1:1:n_max
Myz(:,i) = Myz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl=dlmread ('Mzx.txt', '' ,5,0);
Mzx  = ones(length_fre, n_max);
for i = 1:1:n_max
Mzx(:,i) = Mzx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Mzy.txt', '' ,5,0);
Mzy  = ones(length_fre, n_max);
for i = 1:1:n_max
Mzy(:,i) = Mzy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Mzz.txt', '' ,5,0);
Mzz  = ones(length_fre, n_max);
for i = 1:1:n_max
Mzz(:,i) = Mzz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qxx.txt', '' ,5,0);
Qxx  = ones(length_fre, n_max);
for i = 1:1:n_max
Qxx(:,i) = Qxx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qyy.txt', '' ,5,0);
Qyy  = ones(length_fre, n_max);
for i = 1:1:n_max
Qyy(:,i) = Qyy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qzz.txt', '' ,5,0);
Qzz  = ones(length_fre, n_max);
for i = 1:1:n_max
Qzz(:,i) = Qzz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qxy.txt', '' ,5,0);
Qxy  = ones(length_fre, n_max);
for i = 1:1:n_max
Qxy(:,i) = Qxy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qxz.txt', '' ,5,0);
Qxz  = ones(length_fre, n_max);
for i = 1:1:n_max
Qxz(:,i) = Qxz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qyx.txt', '' ,5,0);
Qyx  = ones(length_fre, n_max);
for i = 1:1:n_max
Qyx(:,i) = Qyx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qyz.txt', '' ,5,0);
Qyz  = ones(length_fre, n_max);
for i = 1:1:n_max
Qyz(:,i) = Qyz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qzx.txt', '' ,5,0);
Qzx  = ones(length_fre, n_max);
for i = 1:1:n_max
Qzx(:,i) = Qzx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Qzy.txt', '' ,5,0);
Qzy  = ones(length_fre, n_max);
for i = 1:1:n_max
Qzy(:,i) = Qzy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oxxx.txt', '' ,5,0);
Oxxx  = ones(length_fre, n_max);
for i = 1:1:n_max
Oxxx(:,i) = Oxxx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oxxy.txt', '' ,5,0);
Oxxy  = ones(length_fre, n_max);
for i = 1:1:n_max
Oxxy(:,i) = Oxxy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oxxz.txt', '' ,5,0);
Oxxz  = ones(length_fre, n_max);
for i = 1:1:n_max
Oxxz(:,i) = Oxxz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oyyx.txt', '' ,5,0);
Oyyx  = ones(length_fre, n_max);
for i = 1:1:n_max
Oyyx(:,i) = Oyyx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oyyy.txt', '' ,5,0);
Oyyy  = ones(length_fre, n_max);
for i = 1:1:n_max
Oyyy(:,i) = Oyyy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oyyz.txt', '' ,5,0);
Oyyz  = ones(length_fre, n_max);
for i = 1:1:n_max
Oyyz(:,i) = Oyyz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Ozzx.txt', '' ,5,0);
Ozzx  = ones(length_fre, n_max);
for i = 1:1:n_max
Ozzx(:,i) = Ozzx(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Ozzy.txt', '' ,5,0);
Ozzy  = ones(length_fre, n_max);
for i = 1:1:n_max
Ozzy(:,i) = Ozzy(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Ozzz.txt', '' ,5,0);
Ozzz  = ones(length_fre, n_max);
for i = 1:1:n_max
Ozzz(:,i) = Ozzz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Oxyz.txt', '' ,5,0);
Oxyz  = ones(length_fre, n_max);
for i = 1:1:n_max
Oxyz(:,i) = Oxyz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Lambdax.txt', '' ,5,0);
Lambdax  = ones(length_fre, n_max);
for i = 1:1:n_max
Lambdax(:,i) = Lambdax(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Lambday.txt', '' ,5,0);
Lambday  = ones(length_fre, n_max);
for i = 1:1:n_max
Lambday(:,i) = Lambday(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('Lambdaz.txt', '' ,5,0);
Lambdaz  = ones(length_fre, n_max);
for i = 1:1:n_max
Lambdaz(:,i) = Lambdaz(:,i) .* epsilon_tbl(i:n_max:end,4);
end

epsilon_tbl = dlmread ('absCS.txt', '' ,5,0); % поглощение в Ваттах
absCS  = ones(length_fre, n_max);
for i = 1:1:n_max
absCS(:,i) = absCS(:,i) .* epsilon_tbl(i:n_max:end,4);
end
absCS = absCS ./ 1.3E-3;    % нормировка поглощения на вектор пойнтинга падающего излучения

% Компоненты экстинкции для волны направленной по Z и поляризованной по X
ExtPx 	= ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* Px);
ExtTx 	= ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* ((1i.*kd./vd).*1.*Tx) );
ExtQxz	= ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (-(1i.*kd./6).*Qxz) );
Extmy 	= ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (my./vd) );
ExtMyz 	= ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* ((1i.*kd./(2.*vd)).*Myz) );
ExtOxzz = ( kd ./ (eps0.*epsd .* abs(E0x).^2) ) .* imag( conj(E0x) .* (-(kd.*kd./6) .* (Ozzx-Lambdax)) );


TxK = ((1i.*kd./vd) .* Tx);
TyK = ((1i.*kd./vd) .* Ty);
TzK = ((1i.*kd./vd) .* Tz);

% Экстинкция как сумма мультипольных составляющий
ExtCS = ExtPx+ExtTx+ExtQxz+Extmy+ExtMyz+ExtOxzz;

% Компоненты полного электрического дипольного момента
Dx = Px - (1i.*k0.*epsd./c).*Tx;
Dy = Py - (1i.*k0.*epsd./c).*Ty;
Dz = Pz - (1i.*k0.*epsd./c).*Tz;

% I incident - вектор пойнтинга падающей волны, плотность потока мощности
Iinc = ((eps0*epsd/mu0)^0.5*E0/2);

% Вклады отдельных мультиполей в рассеяние мощности нормированные на вектор пойнтинга
ScatD = (k0.^4./(12.*pi.*eps0.^2.*vd2.*mu0)) .* (abs(Dx).^2 + abs(Dy).^2 + abs(Dz).^2) ./ Iinc; 	% диполь

Scatm = ((k0.^4.*epsd./(12.*pi*eps0.*vd2)) .* (abs(mx).^2+abs(my).^2+abs(mz).^2)) ./ Iinc;  % магнитный диполь

ScatQ=((k0.^6.*epsd./(1440.*pi.*eps0.^2.*vd2.*mu0))... 									% электрический квадруполь
    .*(abs(Qxx).^2+abs(Qxy).^2+abs(Qxz).^2+abs(Qyx).^2+abs(Qyy).^2+abs(Qyz).^2+...
    abs(Qzx).^2+abs(Qzy).^2+abs(Qzz).^2)) ./ Iinc; 											

ScatM = ((k0.^6.*epsd.^2./(160.*pi*eps0.^1.*vd2)) ...									% магнитный квадруполь
    .*(abs(Mxx).^2+abs(Mxy).^2+abs(Mxz).^2+abs(Myx).^2+abs(Myy).^2+abs(Myz).^2+...
    abs(Mzx).^2+abs(Mzy).^2+abs(Mzz).^2) ) ./ Iinc; 										

% запчасти от октуполья из комсола
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

ScatO = ((k0.^8.*epsd.^2./(3780.*pi.*eps0.^2.*vd2.*mu0)).*(6.*abs(FullOxyz).^2+3.*abs(FullOxxy).^2+...	% электрический октуполь
    3.*abs(FullOxxz).^2+3.*abs(FullOyyx).^2+3.*abs(FullOyyz).^2+3.*abs(FullOzzx).^2+...
    3.*abs(FullOzzy).^2+abs(FullOxxx).^2+abs(FullOyyy).^2+abs(FullOzzz).^2)) ./ Iinc ;  		
                
% Сечение рассеяние как сумма мультипольных компонентов
ScatCS = (ScatD + Scatm + ScatQ + ScatM + ScatO) ; 


lambda_nm = lambda .* norm_length; 	% переводим длину волны из метров в нанометры
H = H .* norm_length; 				% переводим высоту в нанометры
FontSize = 15; 			% размер шрифта в подписях
LineWidth = 2.3; 		% толщина линии на графиках

fig1 = figure (1);
% hold on;
pl1 = @(n) plot (lambda_nm(:,n), abs(Px(:,n))./geomCS, ...
	 			lambda_nm(:,n), abs(TxK(:,n) + Px(:,n))./geomCS, ...
	 			lambda_nm(:,n), abs(TxK(:,n))./geomCS, ...
	 			'LineWidth', LineWidth);
xlab1 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab1 = @() ylabel ('dipole, a.u.','FontSize', FontSize);
leg1 =  @() legend( 'ED','ED+TD','TD');
tit1 = 	@(n) title(strcat('Height h = ', num2str(H(1,n)), ' nm'),'FontSize', FontSize);
slider_toroidal( fig1, pl1, xlab1, ylab1, leg1, tit1, n, n_min, n_max, H);
% hold off;


% Multipoles contributions to Extincntion
fig2 = figure(2);
% hold on;
pl2 = @(n) plot (lambda_nm(:,n), ExtPx(:,n)./geomCS, ...
				lambda_nm(:,n), ExtTx(:,n)./geomCS, ...
				lambda_nm(:,n), ExtQxz(:,n)./geomCS, ...
				lambda_nm(:,n), Extmy(:,n)./geomCS, ...
				lambda_nm(:,n), ExtMyz(:,n)./geomCS, ...
				lambda_nm(:,n), ExtOxzz(:,n)./geomCS, ...
                lambda_nm(:,n), ExtCS(:,n)./geomCS, ...
	 			'LineWidth', LineWidth);
tit2 = @(n) title(strcat('Multipoles contributions to Extincntion, h = ', num2str(H(1,n)), ' nm' ));
xlab2 = @() xlabel ('Wavelength ,nm','FontSize', FontSize);
ylab2 = @() ylabel ('Multipoles contributions, um^2','FontSize', FontSize); % um - микрометры
leg2 = @() legend('Px','Tx','Qxz','my','Myz','Oxzz','Total Ext-on Cross Sect. as Sum');
slider_toroidal( fig2, pl2, xlab2, ylab2, leg2, tit2, n, n_min, n_max, H );
% hold off;


% Все сечения 
fig3 = figure(3);
% hold on;
pl3 = @(n) plot (lambda_nm(:,n), abs(ExtCS(:,n))./geomCS, ...
				lambda_nm(:,n), scat(:,n)./geomCS, ...
				lambda_nm(:,n), absCS(:,n)./geomCS, ...
				lambda_nm(:,n), (scat(:,n)+absCS(:,n))./geomCS, ...
				lambda_nm(:,n), ScatCS(:,n)./geomCS, ...
	 			'LineWidth', LineWidth);
tit3 = @(n) title(strcat('Cross section, h = ', num2str(H(1,n)), ' nm' ));
xlab3 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab3 = @() ylabel ('Cross section, mkm^2','FontSize', FontSize);
leg3 = @() legend('Extincntion cross section','Scattering cross section from COMSOL','Absorption cross section from COMSOL', 'Extinction cross section from COMSOL (Abs+Scat)', 'Scatterning cross section' );
slider_toroidal( fig3, pl3, xlab3, ylab3, leg3, tit3, n, n_min, n_max, H );
% hold off;


% Только рассеяния
fig4 = figure (4);
% hold on;
pl4 = @(n) plot (lambda_nm(:,n), ScatCS(:,n)./geomCS, ...
				lambda_nm(:,n), scat(:,n)./geomCS, ...
	 			'LineWidth', LineWidth);
tit4 = @(n) title(strcat('Cross section, h = ', num2str(H(1,n)),' nm' ));
xlab4 = @() xlabel ('Wavelength, nm','FontSize', FontSize);
ylab4 = @() ylabel ('Cross section, um^2','FontSize', FontSize);
leg4 = @() legend('Scattering cross section', 'Scattering cross section from COMSOL');
slider_toroidal( fig4, pl4, xlab4, ylab4, leg4, tit4, n, n_min, n_max, H );
% hold off;


fig5 = figure (5);
% hold on;
pl5 = @(n) plot (lambda_nm(:,n), ScatD(:,n)./geomCS, ...
				lambda_nm(:,n), Scatm(:,n)./geomCS, ...
				lambda_nm(:,n), ScatQ(:,n)./geomCS, ...
				lambda_nm(:,n), ScatM(:,n)./geomCS, ...
				lambda_nm(:,n), ScatO(:,n)./geomCS, ...
				lambda_nm(:,n), ScatCS(:,n)./geomCS, ...
	 			'LineWidth', LineWidth);
tit5 = @(n) title(strcat('Multipoles Contributions to Scattering, h = ', num2str(H(1,n)), ' nm' ));
xlab5 = @() xlabel ('Wavelenght, nm','FontSize', FontSize);
ylab5 = @() ylabel ('Multipoles Contributions, um^2','FontSize', FontSize);
leg5 = @() legend('scat D ', 'scat m', 'scat Q', 'scat M', 'scat O', 'Sum Scat');
slider_toroidal( fig5, pl5, xlab5, ylab5, leg5, tit5, n, n_min, n_max, H );
% hold off;


% %{
% MaxValue [4, n_max]
% 1 строка - максимальное значение величины
% 2 строка - частота соответствующая этому значению
% 3 строка - длина волны соответствующая этому значению
% 4 строка - значение параметра которому соответствует максимальное значение
MaxTxk = MaxValue(abs(TxK), fre, H, n_max); 

figure(6);
% Зависимость максимального значения Txk от высоты
subplot(2,2,1)
plot(MaxTxk(4,:), MaxTxk(1,:), 'o');

% Зависимость частоты соответствующей максимальному значению Txk от высоты
subplot(2,2,3)
plot(MaxTxk(4,:), MaxTxk(2,:), 'o');

% Зависимость длины волны соответствующей максимальному значению Txk от высоты
subplot(2,2,4)
plot(MaxTxk(4,:), MaxTxk(3,:), 'o'); 
% %}

%{
figure;
% hold on;
plot (lambda.*1e9, abs(ExtPx)./geomCS,'b', lambda.*1e9, abs(ExtTx)./geomCS,'y',lambda.*1e9, abs(ExtQxz)./geomCS,'g',lambda.*1e9, abs(Extmy)./geomCS,'m',lambda.*1e9, abs(ExtMyz)./geomCS,'r',lambda.*1e9, abs(ExtOxzz)./geomCS,'k',...
    lambda.*1e9, ScatD./geomCS,'b.', lambda.*1e9, Scatm./geomCS,'m.', lambda.*1e9, ScatQ./geomCS,'g.', lambda.*1e9, ScatM./geomCS,'r.', lambda.*1e9, ScatO./geomCS,'k.');  
xlabel ('Wavelength, nm','FontSize',15);
ylabel ('Multipoles Contributions, mkm^2','FontSize',15);
legend('Px','Tx','Qxz','my','Myz','Oxzz', 'scat D ', 'scat m', 'scat Q', 'scat M', 'scat O');
% hold off; 
%}


