% Allows to delete calculation for 1 parameter value.
% ver 5.1
% 19.03.2017

%n = 3 % Номер расчёта, который надо удалить

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

button = questdlg('Are you sure do you want to delete calculation?','Adequacy test','Yes','No','Yes');
switch button
	case 'Yes'

		prompt = {'Enter calculation index number to delete:'};
		dlg_title = 'Calculation to delete';
		num_lines = 1;
		defaultans = {'NaN'};
		n = inputdlg(prompt,dlg_title,num_lines,defaultans);
        n = str2num(n{1}(1));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Считываем файлы

		fre     = dlmread ('fre.dat');
		H       = dlmread ('H.dat');
		Px      = dlmread ('Px.dat');
		Py      = dlmread ('Py.dat');
		Pz      = dlmread ('Pz.dat');
		Tx      = dlmread ('Tx.dat');
		Ty      = dlmread ('Ty.dat');
		Tz      = dlmread ('Tz.dat');
		mx      = dlmread ('mx.dat');
		my      = dlmread ('my.dat');
		mz      = dlmread ('mz.dat');
		Mxx     = dlmread ('Mxx.dat');
		Mxy     = dlmread ('Mxy.dat');
		Mxz     = dlmread ('Mxz.dat');
		Myx     = dlmread ('Myx.dat');
		Myy     = dlmread ('Myy.dat');
		Myz     = dlmread ('Myz.dat');
		Mzx     = dlmread ('Mzx.dat');
		Mzy     = dlmread ('Mzy.dat');
		Mzz     = dlmread ('Mzz.dat');
		Qxx     = dlmread ('Qxx.dat');
		Qyy     = dlmread ('Qyy.dat');
		Qzz     = dlmread ('Qzz.dat');
		Qxy     = dlmread ('Qxy.dat');
		Qxz     = dlmread ('Qxz.dat');
		Qyx     = dlmread ('Qyx.dat');
		Qyz     = dlmread ('Qyz.dat');
		Qzx     = dlmread ('Qzx.dat');
		Qzy     = dlmread ('Qzy.dat');
		Oxxx    = dlmread ('Oxxx.dat');
		Oxxy    = dlmread ('Oxxy.dat');
		Oxxz    = dlmread ('Oxxz.dat');
		Oyyx    = dlmread ('Oyyx.dat');
		Oyyy    = dlmread ('Oyyy.dat');
		Oyyz    = dlmread ('Oyyz.dat');
		Ozzx    = dlmread ('Ozzx.dat');
		Ozzy    = dlmread ('Ozzy.dat');
		Ozzz    = dlmread ('Ozzz.dat');
		Oxyz    = dlmread ('Oxyz.dat');
		Lambdax = dlmread ('Lambdax.dat');
		Lambday = dlmread ('Lambday.dat');
		Lambdaz = dlmread ('Lambdaz.dat');
		absCS   = dlmread ('absCS.dat');
		scat    = dlmread ('scat.dat');
		ForScat = dlmread ('ForScat.dat');
		BackScat = dlmread ('BackScat.dat');
		ForScatPoint = dlmread ('ForScatPoint.dat');
		BackScatPoint = dlmread ('BackScatPoint.dat');

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Удаляем ненужную строку

		fre (n,:) = [];
		H (n,:) = [];
		Px (n,:) = []; 		
		Py (n,:) = []; 		
		Pz (n,:) = []; 		
		Tx (n,:) = []; 		
		Ty (n,:) = []; 		
		Tz (n,:) = []; 		
		mx (n,:) = []; 		
		my (n,:) = []; 		
		mz (n,:) = []; 		
		Mxx (n,:) = []; 	
		Mxy (n,:) = []; 	
		Mxz (n,:) = []; 	
		Myx (n,:) = []; 	
		Myy (n,:) = []; 	
		Myz (n,:) = []; 	
		Mzx (n,:) = []; 	
		Mzy (n,:) = []; 	
		Mzz (n,:) = []; 	
		Qxx (n,:) = []; 	
		Qyy (n,:) = []; 	
		Qzz (n,:) = []; 	
		Qxy (n,:) = []; 	
		Qxz (n,:) = []; 	
		Qyx (n,:) = []; 	
		Qyz (n,:) = []; 	
		Qzx (n,:) = []; 	
		Qzy (n,:) = []; 	
		Oxxx (n,:) = []; 	
		Oxxy (n,:) = []; 	
		Oxxz (n,:) = []; 	
		Oyyx (n,:) = []; 	
		Oyyy (n,:) = []; 	
		Oyyz (n,:) = []; 	
		Ozzx (n,:) = []; 	
		Ozzy (n,:) = []; 	
		Ozzz (n,:) = []; 	
		Oxyz (n,:) = []; 	
		Lambdax (n,:) = []; 
		Lambday (n,:) = []; 
		Lambdaz (n,:) = []; 
		scat (n,:) = []; 	
		absCS (n,:) = []; 	
		ForScat (n,:) = [];
		BackScat (n,:) = [];
		ForScatPoint (n,:) = [];
		BackScatPoint (n,:) = [];

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Записываем файлы

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

		buttonYes = questdlg(strcat('Calculation №', num2str(n),' has been deleted'),'Success!','Ok','Ok');

	case 'No'
		buttonNo = questdlg('Operation has been canceled','Cancel','Ok','Ok');
end
