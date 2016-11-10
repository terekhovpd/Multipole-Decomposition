% Data concatenation
% ver 3.0

clc
clear all;
 
n = 2;     	             

fre = cell(n,1);
	for i =1:n
		fre{i} = dlmread (strcat('fre', num2str(i), '.dat'));
	end

% Height, Meters
% Высота в метрах
H = cell(n,1);
	for i =1:n
		H{i} = dlmread (strcat('H', num2str(i), '.dat'));
	end

Px = cell(n,1);
	for i = 1:n
		Px{i} = dlmread (strcat('Px', num2str(i), '.dat'));
	end

Py = cell(n,1);
	for i =1:n
		Py{i} = dlmread (strcat('Py', num2str(i), '.dat'));
	end

Pz = cell(n,1);
	for i =1:n
		Pz{i} = dlmread (strcat('Pz', num2str(i), '.dat'));
	end

Px = cell(n,1);
	for i =1:n
		Px{i} = dlmread (strcat('Px', num2str(i), '.dat'));
	end

Tx = cell(n,1);
	for i =1:n
		Tx{i} = dlmread (strcat('Tx', num2str(i), '.dat'));
	end

Ty = cell(n,1);
	for i =1:n
		Ty{i} = dlmread (strcat('Ty', num2str(i), '.dat'));
	end

Tz = cell(n,1);
	for i =1:n
		Tz{i} = dlmread (strcat('Tz', num2str(i), '.dat'));
	end

scat = cell(n,1);
	for i =1:n
		scat{i} = dlmread (strcat('scat', num2str(i), '.dat'));
	end

mx = cell(n,1);
	for i =1:n
		mx{i} = dlmread (strcat('mx', num2str(i), '.dat'));
	end

my = cell(n,1);
	for i =1:n
		my{i} = dlmread (strcat('my', num2str(i), '.dat'));
	end

mz = cell(n,1);
	for i =1:n
		mz{i} = dlmread (strcat('mz', num2str(i), '.dat'));
	end

Mxx = cell(n,1);
	for i =1:n
		Mxx{i} = dlmread (strcat('Mxx', num2str(i), '.dat'));
	end

Mxy = cell(n,1);
	for i =1:n
		Mxy{i} = dlmread (strcat('Mxy', num2str(i), '.dat'));
	end

Mxz = cell(n,1);
	for i =1:n
		Mxz{i} = dlmread (strcat('Mxz', num2str(i), '.dat'));
	end

Myx= cell(n,1);
	for i =1:n
		Myx{i} = dlmread (strcat('Myx', num2str(i), '.dat'));
	end

Myy = cell(n,1);
	for i =1:n
		Myy{i} = dlmread (strcat('Myy', num2str(i), '.dat'));
	end

Myz = cell(n,1);
	for i =1:n
		Myz{i} = dlmread (strcat('Myz', num2str(i), '.dat'));
	end

Mzx = cell(n,1);
	for i =1:n
		Mzx{i} = dlmread (strcat('Mzx', num2str(i), '.dat'));
	end

Mzy = cell(n,1);
	for i =1:n
		Mzy{i} = dlmread (strcat('Mzy', num2str(i), '.dat'));
	end

Mzz = cell(n,1);
	for i =1:n
		Mzz{i} = dlmread (strcat('Mzz', num2str(i), '.dat'));
	end

Qxx = cell(n,1);
	for i =1:n
		Qxx{i} = dlmread (strcat('Qxx', num2str(i), '.dat'));
	end

Qyy = cell(n,1);
	for i =1:n
		Qyy{i} = dlmread (strcat('Qyy', num2str(i), '.dat'));
	end

Qzz = cell(n,1);
	for i =1:n
		Qzz{i} = dlmread (strcat('Qzz', num2str(i), '.dat'));
	end

Qxy = cell(n,1);
	for i =1:n
		Qxy{i} = dlmread (strcat('Qxy', num2str(i), '.dat'));
	end

Qxz = cell(n,1);
	for i =1:n
		Qxz{i} = dlmread (strcat('Qxz', num2str(i), '.dat'));
	end

Qyx = cell(n,1);
	for i =1:n
		Qyx{i} = dlmread (strcat('Qyx', num2str(i), '.dat'));
	end

Qyz = cell(n,1);
	for i =1:n
		Qyz{i} = dlmread (strcat('Qyz', num2str(i), '.dat'));
	end

Qzx = cell(n,1);
	for i =1:n
		Qzx{i} = dlmread (strcat('Qzx', num2str(i), '.dat'));
	end

Qzy = cell(n,1);
	for i =1:n
		Qzy{i} = dlmread (strcat('Qzy', num2str(i), '.dat'));
	end

Oxxx = cell(n,1);
	for i =1:n
		Oxxx{i} = dlmread (strcat('Oxxx', num2str(i), '.dat'));
	end

Oxxy = cell(n,1);
	for i =1:n
		Oxxy{i} = dlmread (strcat('Oxxy', num2str(i), '.dat'));
	end

Oxxz = cell(n,1);
	for i =1:n
		Oxxz{i} = dlmread (strcat('Oxxz', num2str(i), '.dat'));
	end

Oyyx = cell(n,1);
	for i =1:n
		Oyyx{i} = dlmread (strcat('Oyyx', num2str(i), '.dat'));
	end

Oyyy = cell(n,1);
	for i =1:n
		Oyyy{i} = dlmread (strcat('Oyyy', num2str(i), '.dat'));
	end

Oyyz = cell(n,1);
	for i =1:n
		Oyyz{i} = dlmread (strcat('Oyyz', num2str(i), '.dat'));
	end

Ozzx = cell(n,1);
	for i =1:n
		Ozzx{i} = dlmread (strcat('Ozzx', num2str(i), '.dat'));
	end

Ozzy = cell(n,1);
	for i =1:n
		Ozzy{i} = dlmread (strcat('Ozzy', num2str(i), '.dat'));
	end

Ozzz = cell(n,1);
	for i =1:n
		Ozzz{i} = dlmread (strcat('Ozzz', num2str(i), '.dat'));
	end

Oxyz = cell(n,1);
	for i =1:n
		Oxyz{i} = dlmread (strcat('Oxyz', num2str(i), '.dat'));
	end

Lambdax = cell(n,1);
	for i =1:n
		Lambdax{i} = dlmread (strcat('Lambdax', num2str(i), '.dat'));
	end

Lambday = cell(n,1);
	for i =1:n
		Lambday{i} = dlmread (strcat('Lambday', num2str(i), '.dat'));
	end

Lambdaz = cell(n,1);
	for i =1:n
		Lambdaz{i} = dlmread (strcat('Lambdaz', num2str(i), '.dat'));
	end

absCS = cell(n,1);
	for i =1:n
		absCS{i} = dlmread (strcat('absCS', num2str(i), '.dat'));
	end

ForScat = cell(n,1);
	for i =1:n
		ForScat{i} = dlmread (strcat('ForScat', num2str(i), '.dat'));
	end	

BackScat = cell(n,1);
	for i =1:n
		BackScat{i} = dlmread (strcat('BackScat', num2str(i), '.dat'));
	end		

ForScatPoint = cell(n,1);
	for i =1:n
		ForScatPoint{i} = dlmread (strcat('ForScatPoint', num2str(i), '.dat'));
	end	

BackScatPoint = cell(n,1);
	for i =1:n
		BackScatPoint{i} = dlmread (strcat('BackScatPoint', num2str(i), '.dat'));
	end		

clear i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Сшиваем массивы частот и делаем маску
mask = NaN( size(fre{1},1) + size(fre{2}, 1), size(fre{1},2) + size(fre{2}, 2) );
FRE = NaN( 1, size(fre{1},2) + size(fre{2}, 2) );

n = 1;
i = 1;
j = 1;
while n 
	if fre{1}(1,i) < fre{2}(1,j)
		FRE(n) = fre{1}(1,i);
		mask( 1 : size(fre{1}, 1), n ) = 1;
		mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n ) = 0;
		i = i+1;
		if i > size(fre{1},2)  
			FRE( n+1 : n + size( fre{2}(1, j:end), 2) ) = fre{2}(1, j:end);
			mask( 1 : size(fre{1}, 1), n+1 : n + size( fre{2}(1, j:end), 2) ) = 0;
			mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1),  n+1 : n + size( fre{2}(1, j:end), 2) ) = 1;
			break
		end
	elseif fre{1}(1,i) == fre{2}(1,j)
		FRE(n) = fre{1}(1,i);
		mask( 1 : size(fre{1}, 1), n ) = 1;
		mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n ) = 1;
		i = i+1;
		j = j+1;
		if i > size(fre{1},2)  &&  j > size(fre{2},2)
			break
		elseif i > size(fre{1},2)
			FRE( n+1 : n + size( fre{2}(1, j:end), 2) ) = fre{2}(1, j:end);
			mask( 1 : size(fre{1}, 1), n+1 : n + size( fre{2}(1, j:end), 2) ) = 0;
			mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n+1 : n + size( fre{2}(1, j:end), 2) ) = 1;
			break
		elseif j > size(fre{2},2)
			FRE( n+1 : n + size( fre{1}(1, i:end), 2) ) = fre{1}(1, i:end);
			mask( 1 : size(fre{1}, 1), n+1 : n + size( fre{1}(1, i:end), 2) ) = 1;
			mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n+1 : n + size( fre{1}(1, i:end), 2) ) = 0;
			break
		end
	else
		FRE(n) = fre{2}(1,j);
		mask( 1 : size(fre{1}, 1), n ) = 0;
		mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n ) = 1;
		j = j+1;
		if j > size(fre{2},2)
			FRE( n+1 : n + size( fre{1}(1, i:end), 2) ) = fre{1}(1, i:end);
			mask( 1 : size(fre{1}, 1)  ,  n+1 : n + size( fre{1}(1, i:end), 2) ) = 1;
			mask( size(fre{1}, 1) + 1  :  size(fre{1}, 1) + size(fre{2}, 1), n+1 : n + size( fre{1}(1, i:end), 2) ) = 0;
			break
		end
	end
	n = n+1;
end

i=0;
while i < size(mask, 2) && ~isnan(mask(1,i+1))
	i = i+1;
end
mask = mask(:, 1:i);
FRE = FRE(:, 1:i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Сшиваем все данные
Px 		= Mask( Px, fre, mask );
Py 		= Mask( Py, fre, mask );
Pz 		= Mask( Pz, fre, mask );
Tx 		= Mask( Tx, fre, mask );
Ty 		= Mask( Ty, fre, mask );
Tz 		= Mask( Tz, fre, mask );
mx 		= Mask( mx, fre, mask );
my 		= Mask( my, fre, mask );
mz 		= Mask( mz, fre, mask );
Mxx 	= Mask( Mxx, fre, mask );
Mxy 	= Mask( Mxy, fre, mask );
Mxz 	= Mask( Mxz, fre, mask );
Myx 	= Mask( Myx, fre, mask );
Myy 	= Mask( Myy, fre, mask );
Myz 	= Mask( Myz, fre, mask );
Mzx 	= Mask( Mzx, fre, mask );
Mzy 	= Mask( Mzy, fre, mask );
Mzz 	= Mask( Mzz, fre, mask );
Qxx 	= Mask( Qxx, fre, mask );
Qyy 	= Mask( Qyy, fre, mask );
Qzz 	= Mask( Qzz, fre, mask );
Qxy 	= Mask( Qxy, fre, mask );
Qxz 	= Mask( Qxz, fre, mask );
Qyx 	= Mask( Qyx, fre, mask );
Qyz 	= Mask( Qyz, fre, mask );
Qzx 	= Mask( Qzx, fre, mask );
Qzy 	= Mask( Qzy, fre, mask );
Oxxx 	= Mask( Oxxx, fre, mask );
Oxxy 	= Mask( Oxxy, fre, mask );
Oxxz 	= Mask( Oxxz, fre, mask );
Oyyx 	= Mask( Oyyx, fre, mask );
Oyyy 	= Mask( Oyyy, fre, mask );
Oyyz 	= Mask( Oyyz, fre, mask );
Ozzx 	= Mask( Ozzx, fre, mask );
Ozzy 	= Mask( Ozzy, fre, mask );
Ozzz 	= Mask( Ozzz, fre, mask );
Oxyz 	= Mask( Oxyz, fre, mask );
Lambdax = Mask( Lambdax, fre, mask );
Lambday = Mask( Lambday, fre, mask );
Lambdaz = Mask( Lambdaz, fre, mask );
scat 	= Mask( scat, fre, mask );
absCS 	= Mask( absCS, fre, mask );
ForScat = Mask( ForScat,fre,mask );
BackScat = Mask( BackScat,fre,mask );
ForScatPoint = Mask( ForScatPoint,fre,mask );
BackScatPoint = Mask( BackScatPoint,fre,mask );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Делаем вектор высот
H = [ H{1}(:,1); H{2}(:,1) ];
[H, I] = sort(H); % сортируем вектор высот
% Сортеруем все остальные данные в соответствии сортировкой вектора высот
mask 	= mask(I,:);
Px 		= Px (I,:);	
Py 		= Py (I,:);	
Pz 		= Pz (I,:);	
Tx 		= Tx (I,:);	
Ty 		= Ty (I,:);	
Tz 		= Tz (I,:);	
mx 		= mx (I,:);	
my 		= my (I,:);	
mz 		= mz (I,:);	
Mxx 	= Mxx (I,:);
Mxy 	= Mxy (I,:);
Mxz 	= Mxz (I,:);
Myx 	= Myx (I,:);
Myy 	= Myy (I,:);
Myz 	= Myz (I,:);
Mzx 	= Mzx (I,:);
Mzy 	= Mzy (I,:);
Mzz 	= Mzz (I,:);
Qxx 	= Qxx (I,:);
Qyy 	= Qyy (I,:);
Qzz 	= Qzz (I,:);
Qxy 	= Qxy (I,:);
Qxz 	= Qxz (I,:);
Qyx 	= Qyx (I,:);
Qyz 	= Qyz (I,:);
Qzx 	= Qzx (I,:);
Qzy 	= Qzy (I,:);
Oxxx 	= Oxxx (I,:);
Oxxy 	= Oxxy (I,:);
Oxxz 	= Oxxz (I,:);
Oyyx 	= Oyyx (I,:);
Oyyy 	= Oyyy (I,:);
Oyyz 	= Oyyz (I,:);
Ozzx 	= Ozzx (I,:);
Ozzy 	= Ozzy (I,:);
Ozzz 	= Ozzz (I,:);
Oxyz 	= Oxyz (I,:);
Lambdax = Lambdax (I,:);
Lambday = Lambday (I,:);
Lambdaz = Lambdaz (I,:); 
scat 	= scat (I,:);	
absCS 	= absCS	(I,:);
ForScat = ForScat (I,:);
BackScat = BackScat (I,:);
ForScatPoint = ForScatPoint (I,:);
BackScatPoint = BackScatPoint (I,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Удаляем повторяющиеся строки
i=1;
while i < length(H)
	if H(i) == H(i+1)
		H(i+1) = [];
		mask (i,:) = [];
		Px (i,:) = []; 		
		Py (i,:) = []; 		
		Pz (i,:) = []; 		
		Tx (i,:) = []; 		
		Ty (i,:) = []; 		
		Tz (i,:) = []; 		
		mx (i,:) = []; 		
		my (i,:) = []; 		
		mz (i,:) = []; 		
		Mxx (i,:) = []; 	
		Mxy (i,:) = []; 	
		Mxz (i,:) = []; 	
		Myx (i,:) = []; 	
		Myy (i,:) = []; 	
		Myz (i,:) = []; 	
		Mzx (i,:) = []; 	
		Mzy (i,:) = []; 	
		Mzz (i,:) = []; 	
		Qxx (i,:) = []; 	
		Qyy (i,:) = []; 	
		Qzz (i,:) = []; 	
		Qxy (i,:) = []; 	
		Qxz (i,:) = []; 	
		Qyx (i,:) = []; 	
		Qyz (i,:) = []; 	
		Qzx (i,:) = []; 	
		Qzy (i,:) = []; 	
		Oxxx (i,:) = []; 	
		Oxxy (i,:) = []; 	
		Oxxz (i,:) = []; 	
		Oyyx (i,:) = []; 	
		Oyyy (i,:) = []; 	
		Oyyz (i,:) = []; 	
		Ozzx (i,:) = []; 	
		Ozzy (i,:) = []; 	
		Ozzz (i,:) = []; 	
		Oxyz (i,:) = []; 	
		Lambdax (i,:) = []; 
		Lambday (i,:) = []; 
		Lambdaz (i,:) = []; 
		scat (i,:) = []; 	
		absCS (i,:) = []; 
		ForScat (i,:) = [];
		BackScat (i,:) = [];
		ForScatPoint (i,:) = [];
		BackScatPoint (i,:) = [];
		ForScatPow (i,:) = [];
		BackScatPow (i,:) = [];
		ForScatPointPow (i,:) = [];
		BackScatPointPow  (i,:) = [];	
	else
		i = i+1;
	end
end

fre = repmat( FRE, size(mask, 1), 1 );
H = repmat( H, 1, size(mask, 2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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