%ver 5.1
%19.03.2017
% Internal function for data_export


function norm_length = out_of_comsol

d = dialog('Name','Dimension from COMSOL','Units', 'normalized', 'Position', [0.3 0.4 0.4 0.24]);
norm_length = 1e-9;
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
    'String','What H dimension in the .txt files that have been exported from COMSOL? (Какая размерность в txt-файлах, экспортированных из Comsol-а?)');

waitfor(popup)

function setpopup(source,callbackdata)
    switch(source.Value)
        case 1
            norm_length = 1e-9;
        case 2
            norm_length = 1e-6;
        case 3
            norm_length = 1e-3;
        case 4
            norm_length = 1;
    end
end

end