% 
% ver 4.02

clc
clear all;

n = 1;     	             
            
epsilon_tbl = dlmread ('Px.txt', '' ,5,0); 	% read data from file, delete 5 top strings (header in COMSOL export files) 
                                            % считывет данные из файла, обрезает 5 верхних строк
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency array creating
FRE = epsilon_tbl(:, 1)';	% read all elements from 1 column   
                            % считывает из 1 колонки все элементы
i = 1;
while FRE(i) == FRE(i+1)
    i=i+1;
end
n_max = i;  % number of parametric sweep steps, maximum graphic slider parameter 
            % число шагов в parametric sweep, максимальное значение слайдера
FRE = unique(FRE);  % leave in vector FRE only unique elements
                    % оставляет в векторе FRE только уникальные элементы
length_fre = size(FRE,2);	% number of frequency steps, frequency vector length  
                            % число шагов по частоте, длина вектора частот                                                                 
fre = ones(n_max, length_fre);
for i = 1:1:n_max
	fre(i,:) = fre(i,:) .* FRE; 	% create the matrix, containing n_max column-vectors fre, from one column-vector fre (for matrix dimension match)   
                                    % делаем из вектор-столбца fre матрицу, состоящую из n_max вектор-столбцов fre
end
clear FRE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Height, Meters
% Высота в метрах
H = ones(n_max, length_fre);
for i = 1:1:length_fre
H(:,i) = H(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 2);
end


Px = ones(n_max, length_fre);
for i = 1:1:length_fre
Px(:,i) = Px(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Py.txt', '' ,5,0);
Py = ones(n_max, length_fre);
for i = 1:1:length_fre
Py(:,i) = Py(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Pz.txt', '' ,5,0);
Pz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Pz(:,i) = Pz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Tx.txt', '' ,5,0);
Tx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Tx(:,i) = Tx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Ty.txt', '' ,5,0);
Ty  = ones(n_max, length_fre);
for i = 1:1:length_fre
Ty(:,i) = Ty(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Tz.txt', '' ,5,0);
Tz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Tz(:,i) = Tz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

% Scattering cross-section normalization by Poynting vector of incident radiation
% нормировка сечения рассеяния на вектор пойнтинга падающего излучения
epsilon_tbl=dlmread ('scat.txt', '' ,5,0);
scat  = ones(n_max, length_fre);
for i = 1:1:length_fre
scat(:,i) = scat(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4) ./ 1.3E-3;	% Scattering cross-section, Watts 
                                                                                % сечение рассеяния в Ваттах
end                 % cross-section normalization by Poynting vector of incident radiation  
                    % нормировка сечения на вектор пойнтинга падающего излучения

epsilon_tbl=dlmread ('mx.txt', '' ,5,0);
mx  = ones(n_max, length_fre);
for i = 1:1:length_fre
mx(:,i) = mx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('my.txt', '' ,5,0);
my  = ones(n_max, length_fre);
for i = 1:1:length_fre
my(:,i) = my(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('mz.txt', '' ,5,0);
mz  = ones(n_max, length_fre);
for i = 1:1:length_fre
mz(:,i) = mz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Mxx.txt', '' ,5,0);
Mxx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mxx(:,i) = Mxx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Mxy.txt', '' ,5,0);
Mxy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mxy(:,i) = Mxy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Mxz.txt', '' ,5,0);
Mxz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mxz(:,i) = Mxz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Myx.txt', '' ,5,0);
Myx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Myx(:,i) = Myx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Myy.txt', '' ,5,0);
Myy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Myy(:,i) = Myy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Myz.txt', '' ,5,0);
Myz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Myz(:,i) = Myz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl=dlmread ('Mzx.txt', '' ,5,0);
Mzx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mzx(:,i) = Mzx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Mzy.txt', '' ,5,0);
Mzy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mzy(:,i) = Mzy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Mzz.txt', '' ,5,0);
Mzz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Mzz(:,i) = Mzz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qxx.txt', '' ,5,0);
Qxx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qxx(:,i) = Qxx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qyy.txt', '' ,5,0);
Qyy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qyy(:,i) = Qyy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qzz.txt', '' ,5,0);
Qzz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qzz(:,i) = Qzz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qxy.txt', '' ,5,0);
Qxy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qxy(:,i) = Qxy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qxz.txt', '' ,5,0);
Qxz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qxz(:,i) = Qxz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qyx.txt', '' ,5,0);
Qyx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qyx(:,i) = Qyx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qyz.txt', '' ,5,0);
Qyz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qyz(:,i) = Qyz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qzx.txt', '' ,5,0);
Qzx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qzx(:,i) = Qzx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Qzy.txt', '' ,5,0);
Qzy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Qzy(:,i) = Qzy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oxxx.txt', '' ,5,0);
Oxxx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oxxx(:,i) = Oxxx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oxxy.txt', '' ,5,0);
Oxxy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oxxy(:,i) = Oxxy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oxxz.txt', '' ,5,0);
Oxxz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oxxz(:,i) = Oxxz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oyyx.txt', '' ,5,0);
Oyyx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oyyx(:,i) = Oyyx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oyyy.txt', '' ,5,0);
Oyyy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oyyy(:,i) = Oyyy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oyyz.txt', '' ,5,0);
Oyyz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oyyz(:,i) = Oyyz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Ozzx.txt', '' ,5,0);
Ozzx  = ones(n_max, length_fre);
for i = 1:1:length_fre
Ozzx(:,i) = Ozzx(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Ozzy.txt', '' ,5,0);
Ozzy  = ones(n_max, length_fre);
for i = 1:1:length_fre
Ozzy(:,i) = Ozzy(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Ozzz.txt', '' ,5,0);
Ozzz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Ozzz(:,i) = Ozzz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Oxyz.txt', '' ,5,0);
Oxyz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Oxyz(:,i) = Oxyz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end


epsilon_tbl = dlmread ('Lambdax.txt', '' ,5,0);
Lambdax  = ones(n_max, length_fre);
for i = 1:1:length_fre
Lambdax(:,i) = Lambdax(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Lambday.txt', '' ,5,0);
Lambday  = ones(n_max, length_fre);
for i = 1:1:length_fre
Lambday(:,i) = Lambday(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('Lambdaz.txt', '' ,5,0);
Lambdaz  = ones(n_max, length_fre);
for i = 1:1:length_fre
Lambdaz(:,i) = Lambdaz(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end

epsilon_tbl = dlmread ('absCS.txt', '' ,5,0); % Absorption, Watts % поглощение в Ваттах
absCS  = ones(n_max, length_fre);
for i = 1:1:length_fre
absCS(:,i) = absCS(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end
absCS = absCS ./ 1.3E-3;    % absorption normalization by Poynting vector of incident radiation   % нормировка поглощения на вектор пойнтинга падающего излучения


epsilon_tbl = dlmread ('ForScat.txt', '' ,5,0); % Forward Scattering (integrated by half space)
ForScat  = ones(n_max, length_fre);
for i = 1:1:length_fre
ForScat(:,i) = ForScat(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end
ForScat = ForScat ./ 1.3E-3;

epsilon_tbl = dlmread ('BackScat.txt', '' ,5,0); % Backward Scattering (integrated by half space)
BackScat  = ones(n_max, length_fre);
for i = 1:1:length_fre
BackScat(:,i) = BackScat(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end
BackScat = BackScat ./ 1.3E-3;

epsilon_tbl = dlmread ('ForScatPoint.txt', '' ,5,0); % Forward Scattering (point)
ForScatPoint  = ones(n_max, length_fre);
for i = 1:1:length_fre
ForScatPoint(:,i) = ForScatPoint(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end
ForScatPoint = ForScatPoint ./ 1.3E-3;

epsilon_tbl = dlmread ('BackScatPoint.txt', '' ,5,0); % Backward Scattering (point)
BackScatPoint  = ones(n_max, length_fre);
for i = 1:1:length_fre
BackScatPoint(:,i) = BackScatPoint(:,i) .* epsilon_tbl(1+(i-1)*n_max : 1 : i*n_max, 4);
end
BackScatPoint = BackScatPoint ./ 1.3E-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n == 0
        dlmwrite ('fre.dat', fre, 'delimiter', '\t');
        dlmwrite ('H.dat', H, 'delimiter', '\t');
        dlmwrite ('Px.dat', Px, 'delimiter', '\t');
        dlmwrite ('Py.dat', Py, 'delimiter', '\t');
        dlmwrite ('Pz.dat', Pz, 'delimiter', '\t');
        dlmwrite ('Tx.dat', Tx, 'delimiter', '\t');
        dlmwrite ('Ty.dat', Ty, 'delimiter', '\t');
        dlmwrite ('Tz.dat', Tz, 'delimiter', '\t');
        dlmwrite ('mx.dat', mx, 'delimiter', '\t');
        dlmwrite ('my.dat', my, 'delimiter', '\t');
        dlmwrite ('mz.dat', mz, 'delimiter', '\t');
        dlmwrite ('Mxx.dat', Mxx, 'delimiter', '\t');
        dlmwrite ('Mxy.dat', Mxy, 'delimiter', '\t');
        dlmwrite ('Mxz.dat', Mxz, 'delimiter', '\t');
        dlmwrite ('Myx.dat', Myx, 'delimiter', '\t');
        dlmwrite ('Myy.dat', Myy, 'delimiter', '\t');
        dlmwrite ('Myz.dat', Myz, 'delimiter', '\t');
        dlmwrite ('Mzx.dat', Mzx, 'delimiter', '\t');
        dlmwrite ('Mzy.dat', Mzy, 'delimiter', '\t');
        dlmwrite ('Mzz.dat', Mzz, 'delimiter', '\t');
        dlmwrite ('Qxx.dat', Qxx, 'delimiter', '\t');
        dlmwrite ('Qyy.dat', Qyy, 'delimiter', '\t');
        dlmwrite ('Qzz.dat', Qzz, 'delimiter', '\t');
        dlmwrite ('Qxy.dat', Qxy, 'delimiter', '\t');
        dlmwrite ('Qxz.dat', Qxz, 'delimiter', '\t');
        dlmwrite ('Qyx.dat', Qyx, 'delimiter', '\t');
        dlmwrite ('Qyz.dat', Qyz, 'delimiter', '\t');
        dlmwrite ('Qzx.dat', Qzx, 'delimiter', '\t');
        dlmwrite ('Qzy.dat', Qzy, 'delimiter', '\t');
        dlmwrite ('Oxxx.dat', Oxxx, 'delimiter', '\t');
        dlmwrite ('Oxxy.dat', Oxxy, 'delimiter', '\t');
        dlmwrite ('Oxxz.dat', Oxxz, 'delimiter', '\t');
        dlmwrite ('Oyyx.dat', Oyyx, 'delimiter', '\t');
        dlmwrite ('Oyyy.dat', Oyyy, 'delimiter', '\t');
        dlmwrite ('Oyyz.dat', Oyyz, 'delimiter', '\t');
        dlmwrite ('Ozzx.dat', Ozzx, 'delimiter', '\t');
        dlmwrite ('Ozzy.dat', Ozzy, 'delimiter', '\t');
        dlmwrite ('Ozzz.dat', Ozzz, 'delimiter', '\t');
        dlmwrite ('Oxyz.dat', Oxyz, 'delimiter', '\t');
        dlmwrite ('Lambdax.dat', Lambdax, 'delimiter', '\t');
        dlmwrite ('Lambday.dat', Lambday, 'delimiter', '\t');
        dlmwrite ('Lambdaz.dat', Lambdaz, 'delimiter', '\t');
        dlmwrite ('absCS.dat', absCS, 'delimiter', '\t');
        dlmwrite ('scat.dat', scat, 'delimiter', '\t');
        dlmwrite ('ForScat.dat', ForScat, 'delimiter', '\t');
        dlmwrite ('BackScat.dat', BackScat, 'delimiter', '\t');
        dlmwrite ('ForScatPoint.dat', ForScatPoint, 'delimiter', '\t');
        dlmwrite ('BackScatPoint.dat', BackScatPoint, 'delimiter', '\t');
else
        dlmwrite (strcat('fre', num2str(n), '.dat'), fre, 'delimiter', '\t');
        dlmwrite (strcat('H', num2str(n), '.dat'), H, 'delimiter', '\t');
        dlmwrite (strcat('Px', num2str(n), '.dat'), Px, 'delimiter', '\t');
        dlmwrite (strcat('Py', num2str(n), '.dat'), Py, 'delimiter', '\t');
        dlmwrite (strcat('Pz', num2str(n), '.dat'), Pz, 'delimiter', '\t');
        dlmwrite (strcat('Tx', num2str(n), '.dat'), Tx, 'delimiter', '\t');
        dlmwrite (strcat('Ty', num2str(n), '.dat'), Ty, 'delimiter', '\t');
        dlmwrite (strcat('Tz', num2str(n), '.dat'), Tz, 'delimiter', '\t');
        dlmwrite (strcat('mx', num2str(n), '.dat'), mx, 'delimiter', '\t');
        dlmwrite (strcat('my', num2str(n), '.dat'), my, 'delimiter', '\t');
        dlmwrite (strcat('mz', num2str(n), '.dat'), mz, 'delimiter', '\t');
        dlmwrite (strcat('Mxx', num2str(n), '.dat'), Mxx, 'delimiter', '\t');
        dlmwrite (strcat('Mxy', num2str(n), '.dat'), Mxy, 'delimiter', '\t');
        dlmwrite (strcat('Mxz', num2str(n), '.dat'), Mxz, 'delimiter', '\t');
        dlmwrite (strcat('Myx', num2str(n), '.dat'), Myx, 'delimiter', '\t');
        dlmwrite (strcat('Myy', num2str(n), '.dat'), Myy, 'delimiter', '\t');
        dlmwrite (strcat('Myz', num2str(n), '.dat'), Myz, 'delimiter', '\t');
        dlmwrite (strcat('Mzx', num2str(n), '.dat'), Mzx, 'delimiter', '\t');
        dlmwrite (strcat('Mzy', num2str(n), '.dat'), Mzy, 'delimiter', '\t');
        dlmwrite (strcat('Mzz', num2str(n), '.dat'), Mzz, 'delimiter', '\t');
        dlmwrite (strcat('Qxx', num2str(n), '.dat'), Qxx, 'delimiter', '\t');
        dlmwrite (strcat('Qyy', num2str(n), '.dat'), Qyy, 'delimiter', '\t');
        dlmwrite (strcat('Qzz', num2str(n), '.dat'), Qzz, 'delimiter', '\t');
        dlmwrite (strcat('Qxy', num2str(n), '.dat'), Qxy, 'delimiter', '\t');
        dlmwrite (strcat('Qxz', num2str(n), '.dat'), Qxz, 'delimiter', '\t');
        dlmwrite (strcat('Qyx', num2str(n), '.dat'), Qyx, 'delimiter', '\t');
        dlmwrite (strcat('Qyz', num2str(n), '.dat'), Qyz, 'delimiter', '\t');
        dlmwrite (strcat('Qzx', num2str(n), '.dat'), Qzx, 'delimiter', '\t');
        dlmwrite (strcat('Qzy', num2str(n), '.dat'), Qzy, 'delimiter', '\t');
        dlmwrite (strcat('Oxxx', num2str(n), '.dat'), Oxxx, 'delimiter', '\t');
        dlmwrite (strcat('Oxxy', num2str(n), '.dat'), Oxxy, 'delimiter', '\t');
        dlmwrite (strcat('Oxxz', num2str(n), '.dat'), Oxxz, 'delimiter', '\t');
        dlmwrite (strcat('Oyyx', num2str(n), '.dat'), Oyyx, 'delimiter', '\t');
        dlmwrite (strcat('Oyyy', num2str(n), '.dat'), Oyyy, 'delimiter', '\t');
        dlmwrite (strcat('Oyyz', num2str(n), '.dat'), Oyyz, 'delimiter', '\t');
        dlmwrite (strcat('Ozzx', num2str(n), '.dat'), Ozzx, 'delimiter', '\t');
        dlmwrite (strcat('Ozzy', num2str(n), '.dat'), Ozzy, 'delimiter', '\t');
        dlmwrite (strcat('Ozzz', num2str(n), '.dat'), Ozzz, 'delimiter', '\t');
        dlmwrite (strcat('Oxyz', num2str(n), '.dat'), Oxyz, 'delimiter', '\t');
        dlmwrite (strcat('Lambdax', num2str(n), '.dat'), Lambdax, 'delimiter', '\t');
        dlmwrite (strcat('Lambday', num2str(n), '.dat'), Lambday, 'delimiter', '\t');
        dlmwrite (strcat('Lambdaz', num2str(n), '.dat'), Lambdaz, 'delimiter', '\t');
        dlmwrite (strcat('absCS', num2str(n), '.dat'), absCS, 'delimiter', '\t');
        dlmwrite (strcat('scat', num2str(n), '.dat'), scat, 'delimiter', '\t');
        dlmwrite (strcat('ForScat', num2str(n), '.dat'), ForScat, 'delimiter', '\t');
        dlmwrite (strcat('BackScat', num2str(n), '.dat'), BackScat, 'delimiter', '\t');
        dlmwrite (strcat('ForScatPoint', num2str(n), '.dat'), ForScatPoint, 'delimiter', '\t');
        dlmwrite (strcat('BackScatPoint', num2str(n), '.dat'), BackScatPoint, 'delimiter', '\t');        
end