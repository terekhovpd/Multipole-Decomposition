function slider_toroidal( fig, pl, xlab, ylab, leg, tit, n, n_min, n_max, h )
% 1.4
% Слайдер для тороидальной декомпозиции
% fig - дескриптор фигуры
% pl - функция, которая строит график
% xlab - функция, которая задаёт подпись к оси X
% ylab - функция, которая задаёт подпись к оси Y
% n - начальное значение слайдера
% n_min - минимальное значение слайдера
% n_max - максимальное значение слайдера
% h - параметр по которому варьируем 

% set(0,'DefaultAxesFontSize', 15,'DefaultAxesFontName','Arial');
% set(0,'DefaultTextFontSize', 15,'DefaultTextFontName','Arial'); 

set(fig, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

bgcolor = fig.Color;
% bgcolor = 'white';


% Задаём расположение оси X в окне figure
axes('Parent',fig,'position',[0.13 0.22  0.77 0.7]); % в относительных единицах
pl(n);	 % график функции
xlab();	 % ось X
ylab();	 % ось Y
leg();	 % легенда
tit(n);	 % заголовок
grid on;

% Параметры слайдера
sld_width = 0.6; 				% ширина слайдера
sld_heigh = 0.045; 				% высота слайдера
sld_left = (1 - sld_width)/2; 	% расстояние от левого нижнего края окна до левого нижнего края слайдера
sld_bottom = 0.08;				% расстояние от нижнего края окна до нижнего края слайдера

% значение параметра
% txt2 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
%                 'Position',[0.1, 0.05 - sld_heigh/2 , 2*sld_heigh, 0.7*sld_heigh], ...
%                 'FontUnits','normalized', 'FontSize', 0.85, ...
%                 'String', h(1,n) ,'BackgroundColor','white');

% индикатор под слайдером
txt1 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
                'Position',[0.5 - (0.8*sld_heigh)/2, 0.05 - sld_heigh/2 , 0.8*sld_heigh, 0.7*sld_heigh], ...
                'FontUnits','normalized', 'FontSize', 0.85, ...
                'String', n,'BackgroundColor', bgcolor);
% слайдер
sld = uicontrol('Parent', fig, 'Style','slider', 'Units','normalized', ...
              'Position', [sld_left, sld_bottom, sld_width, sld_heigh], ...
              'value', n, 'min', n_min, 'max', n_max, ...
              'SliderStep', [1/n_max 2*1/n_max], ...
              'Callback', @setvalue);


% левое окошко около слайдера
uicontrol('Parent',fig,'Style','text','Units','normalized', ...
                'Position',[sld_left - 1.5*sld_heigh, sld_bottom, 1.5*sld_heigh, sld_heigh], ...
                'FontUnits','normalized', 'FontSize', 0.7, ...
                'String', h (1,n_min),'BackgroundColor', bgcolor);

% правое окошко около слайдера
uicontrol('Parent',fig,'Style','text','Units','normalized', ...
                'Position',[sld_left + sld_width, sld_bottom, 1.5*sld_heigh , sld_heigh], ...
                'FontUnits','normalized', 'FontSize', 0.7, ...
                'String', h (1, n_max),'BackgroundColor', bgcolor);

	function setvalue (source, dataevent)
	    n = fix(source.Value); 	% сохраняем в n целую часть числа source.Value
	    if source.Value >= n + 0.5
	    	source.Value = n + 1;
	    else
	    	source.Value = n;
	    end
	    pl(source.Value);
 	    xlab();
 	    ylab();
 	    leg();
        tit(source.Value);
 %       txt1.String = h (1, source.Value)
        txt1.String = source.Value;
        grid on;
    end

end

