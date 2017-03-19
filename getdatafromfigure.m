% getdatafromfigure.m is a program that allows you to select a figure and a plot
% or a subplot from which to extract the x and y data.
%
% SYNTAX: [x,y]=getdatafromfigure;
%
% The function returns the data to the stack and can be saved as a .mat file as
% well.
%
% DBE 06/27/00

function [x,y]=getdatafromfigure

disp('Click on the plot or subplot from which to extract data then hit any key...');
pause;
child_handles=get(gca,'children');
num_child=size(child_handles,1);
choice=num_child; % Data is always in the last handle?!?!?!? Maybe assume that it is the longest data of any handle

x=get(child_handles(choice),'xdata');
y=get(child_handles(choice),'ydata');

save_file=questdlg('Do you want to save the stolen data?');
switch save_file
  case 'Yes'
    [filename,pathname]=uiputfile('*.mat','Save Stolen Data As: ');
    if filename
      fname=[pathname filename];
      save (fname,'x','y');
      disp(['Saved x and y data to ', fname]);
    else
      disp('Data not saved (1).');
    end
  otherwise
      disp('Data not saved (2).');
end

plot_data=questdlg('Do you want to plot the data?');
switch plot_data
  case 'Yes'
    figure (3) ; hold on; plot(x,y, 'r--','LineWidth', 2.3 );  %,'MarkerSize',10d
    %title(strcat('MQ contribution to scattering, nm^2' ),'FontSize', 16);
    xlabel ('Wavelength, nm','FontSize', 16);
    ylabel ('Dipole moment magnitude, a.u.','FontSize', 16);
  otherwise
    disp('Data not plotted');
    
    %{
    plot_data=questdlg('Do you want to plot the data?');
switch plot_data
  case 'Yes'
    figure (10) ; hold on; plot(x,y, 'LineWidth', 2.3 );  %,'MarkerSize',10d
    title(strcat('MQ contribution to scattering, nm^2' ),'FontSize', 16);
    xlabel ('Wavelength, nm','FontSize', 16);
    ylabel ('Multipole Contribution, nm^2','FontSize', 16);
  otherwise
    disp('Data not plotted');
    %}
    
end

% Display the stolen data

return;