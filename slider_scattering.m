function slider_scattering( fig, pl, pol1, pol2, bound, ss, n, n_min, n_max, f, f_min, f_max, fre, h, norm_length)

set(fig, 'Units', 'normalized', 'OuterPosition', [0.01 0.045 0.98 0.95]);

bgcolor = fig.Color;

% Задаём расположение оси X в окне figure
% axes('Parent',fig,'position',[0.5 0.22  0.46 0.70]); % в относительных единицах
subplot(10,10,[6 90])
pl(n,f);	 % график функции
title(ss);
bound(n,f);
axis([-max(bound(n,f)) max(bound(n,f)) -max(bound(n,f)) max(bound(n,f)) -max(bound(n,f)) max(bound(n,f))]);
grid on;

subplot(10,10,[1 35])
pol1(n,f)
title('Cross-section in x-z surface')

subplot(10,10,[51 85])
pol2(n,f)
title('Cross-section in y-z surface')

lambda = 2.99792458e8./fre.*norm_length;

% Параметры 1 слайдера
sld_width1 = 0.6; 				% ширина слайдера
sld_heigh1 = 0.02; 				% высота слайдера
sld_left1 = (1 - sld_width1)/2; 	% расстояние от левого нижнего края окна до левого нижнего края слайдера
sld_bottom1 = 0.04;				% расстояние от нижнего края окна до нижнего края слайдера

% Параметры 2 слайдера
sld_width2 = 0.6; 				% ширина слайдера
sld_heigh2 = 0.02; 				% высота слайдера
sld_left2 = (1 - sld_width2)/2; 	% расстояние от левого нижнего края окна до левого нижнего края слайдера
sld_bottom2 = 0.12;				% расстояние от нижнего края окна до нижнего края слайдера

% слайдер
sld = uicontrol('Parent', fig, 'Style','slider', 'Units','normalized', ...
              'Position', [sld_left1, sld_bottom1, sld_width1, sld_heigh1], ...
              'value', n, 'min', n_min, 'max', n_max, ...
              'SliderStep', [1/n_max 2*1/n_max], ...
              'Callback', @setvalueV1);
          
sld2 = uicontrol('Parent', fig, 'Style','slider', 'Units','normalized', ...
              'Position', [sld_left2, sld_bottom2, sld_width2, sld_heigh2], ...
              'value', f, 'min', f_min, 'max', f_max, ...
              'SliderStep', [1/f_max 2*1/f_max], ...
              'Callback', @setvalueV2);
          
txt1 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
                'Position',[0.5 - (0.8*sld_heigh1)/2, 0.03 - sld_heigh1/2 , 1.0*sld_heigh1, 1.0*sld_heigh1], ...
                'FontUnits','normalized', 'FontSize', 0.85, ...
                'String', h (n,1),'BackgroundColor', bgcolor);
            
txt2 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
                'Position',[0.5 - (0.8*sld_heigh2)/2, 0.11 - sld_heigh2/2 , 1.5*sld_heigh2, 1.0*sld_heigh2], ...
                'FontUnits','normalized', 'FontSize', 0.85, ...
                'String', lambda (1,f) ,'BackgroundColor', bgcolor);

%txt3 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
%                'Position',[0.01, 0.01, 0.1, 0.1], ...
%                'FontUnits','normalized', 'FontSize', 0.85, ...
%                'String', fre (1,f) ,'BackgroundColor', bgcolor);
%

%
%txt4 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
%                'Position',[0.01, 0.11, 0.1, 0.2], ...
%                'FontUnits','normalized', 'FontSize', 0.85, ...
%                'String', lambda (1,f) ,'BackgroundColor', bgcolor);
%
%txt5 = uicontrol('Parent',fig,'Style','text','Units','normalized',...
%                'Position',[0.01, 0.21, 0.1, 0.3], ...
%                'FontUnits','normalized', 'FontSize', 0.85, ...
%                'String',  h(n,1) ,'BackgroundColor', bgcolor);
          
	function setvalueV1 (source, dataevent)
	    n = fix(source.Value); 	% сохраняем в n целую часть числа source.Value
	    if source.Value >= n + 0.5
	    	source.Value = n + 1;
	    else
	    	source.Value = n;
        end
        subplot(10,10,[6 90])
        pl(source.Value, f);
        title(ss);
        bound(source.Value,f);
        axis([-max(bound(source.Value,f)) max(bound(source.Value,f)) -max(bound(source.Value,f)) max(bound(source.Value,f)) -max(bound(source.Value,f)) max(bound(source.Value,f))]);
        grid on;
        subplot(10,10,[1 35])
        pol1(source.Value,f)
        title('Cross-section in x-z surface')
        subplot(10,10,[51 85])
        pol2(source.Value,f)
        title('Cross-section in y-z surface')
        txt1.String = h(source.Value,1);
%        txt5.String = h(source.Value,1);
    end

	function setvalueV2 (source, dataevent)
        f = fix(source.Value); 
            if source.Value >= f + 0.5
	    	source.Value = f + 1;
            else
	    	source.Value = f;   
            end
        subplot(10,10,[6 90])
        pl(n, source.Value);
        title(ss);
        bound(n,source.Value);
        axis([-max(bound(n,source.Value)) max(bound(n,source.Value)) -max(bound(n,source.Value)) max(bound(n,source.Value)) -max(bound(n,source.Value)) max(bound(n,source.Value))]);
        grid on;
        subplot(10,10,[1 35])
        pol1(n,source.Value)
        title('Cross-section in x-z surface')
        subplot(10,10,[51 85])
        pol2(n,source.Value)
        title('Cross-section in y-z surface')
        txt2.String = lambda(1,source.Value);
%        txt3.String = fre(1,source.Value);
%        txt4.String = lambda(1,source.Value);
        end

end
