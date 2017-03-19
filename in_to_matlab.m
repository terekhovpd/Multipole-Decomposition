%ver 5.1
%19.03.2017
% Internal function for data_export


function dim = in_to_matlab

d = dialog('Name','Dimension on graphics','Units', 'normalized', 'Position', [0.3 0.4 0.4 0.24]);
dim = cell(2,1);
dim= [{'nm'},{1e9}];
% Create pop-up menu

popup = uicontrol('Parent',d,...
       'Style', 'popup',...
       'String', {'nanometers','micrometers','millimeters','meters'},...
       'Units', 'normalized', ...
       'Position', [0.35 0 0.3 0.5],...
       'Callback', @setpopup);    

% Create push button
btn = uicontrol('Parent',d,...
    'Style', 'pushbutton', 'String', 'Ok',...
    'Units', 'normalized', ...
    'Position', [0.45 0.1 0.1 0.15],...
    'Callback', 'delete(gcf)');       

			
% Add a text uicontrol to label the slider.
txt = uicontrol('Parent',d,...
    'Style','text',...
    'Units', 'normalized', ...
    'Position',[0.15 0.6 0.7 0.3],...
    'String','What H dimension would you like to use on graphics? (Какую размерность вы хотите использовать на графиках?)');

waitfor(popup)

function setpopup(source,callbackdata)
    switch(source.Value)
        case 1
            dim = [{'nm'},{1e9}];
        case 2
            dim = [{'um'},{1e6}];
        case 3
            dim = [{'mm'},{1e3}];
        case 4
            dim = [{'m'},{1}];
    end
end

end