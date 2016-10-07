function [ XX ] = Mask( Px, fre, mask )
% Функция для data concatenation
% применяет маску к файлам
% ver 3.0

XX1 = NaN( size(mask,2), size(fre{1},1) );
%{
[row, col] = find( mask( 1 : size(fre{1},1), :) );
for i = 1 : numel(Px{1})
	XX(row(i), col(i)) = Px{1}(i)';
end
%}
n = find( mask( 1 : size(fre{1},1), :)' );
XX1(n) = Px{1}';
XX1 = XX1';

XX2 = NaN( size(mask,2), size(fre{2},1) );
n = find( mask( size(fre{1},1) +1 : size(fre{1},1) + size(fre{2},1), :)' );
XX2(n) = Px{2}';
XX2 = XX2';

XX = [XX1; XX2];

end

