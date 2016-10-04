function [ MIN ] = LocalMin( Px, n_max )
% Ищет локальные минимумы для каждого параметра (каждой высоты)
% ver 1.0 

MIN = cell(n_max,1);
for j = 1:n_max
	for i = 2 : 1 : size(Px,2)-1
		[M,I] = min([Px(j, i-1), Px(j, i), Px(j, i+1)]);
		if I == 2
			MIN{j} = cat(1, MIN{j}, [i, M]);
		end
	end
end

end

