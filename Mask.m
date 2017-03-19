function [ XX ] = Mask( Px, fre, mask )
% Applay mask to files
% ver 5.1
% 17.03.2017

XX1 = NaN( size(mask,2), size(fre{1},1) );

n = find( mask( 1 : size(fre{1},1), :)' );
XX1(n) = Px{1}';
XX1 = XX1';

XX2 = NaN( size(mask,2), size(fre{2},1) );
n = find( mask( size(fre{1},1) +1 : size(fre{1},1) + size(fre{2},1), :)' );
XX2(n) = Px{2}';
XX2 = XX2';

XX = [XX1; XX2];

end
