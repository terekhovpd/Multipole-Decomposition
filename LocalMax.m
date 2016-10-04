function [ MAX ] = LocalMax( Px, n_max )
% Ищет локальные максимумы для каждого параметра (каждой высоты)
% ver 1.0 

MAX = cell(n_max,1);
for j = 1:n_max
	for i = 2 : 1 : size(Px,2)-1
		[M,I] = max([Px(j, i-1), Px(j, i), Px(j, i+1)]);
		if I == 2
			MAX{j} = cat(1, MAX{j}, [i, M]);
		end
	end
end

end
