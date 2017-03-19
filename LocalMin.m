function [ MIN ] = LocalMin( Px, n_max )
% »щет локальные минимумы дл€ каждого параметра (каждой высоты)
% ver 3.0
% ~2016

temp = cell(n_max, 1); % временный массив €чеек
for i =1:n_max
	temp{i} = Px(i, ~isnan(Px(i,:)) ); % записываем все элементы которые не €вл€ютс€ NaN
end

MIN = cell(n_max,1);
for j = 1:n_max
    for i = 2 : 1 : size(temp{j},2)-1
        [M,I] = min([ temp{j}(i-1), temp{j}(i), temp{j}(i+1) ]);
        if I == 2
            MIN{j} = cat(1, MIN{j}, [find( Px(j,:) == M ), M]);
        end
    end
end

for i = 1:n_max
  if isempty(MIN{i})
      MIN{i} = [1, NaN];
  end
end

end

