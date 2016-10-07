function [ MAX ] = LocalMax( Px, n_max )
% »щет локальные максимумы дл€ каждого параметра (каждой высоты)
% ver 3.0

temp = cell(n_max, 1); % временный массив €чеек
for i =1:n_max
	temp{i} = Px(i, ~isnan(Px(i,:)) ); % записываем все элементы которые не €вл€ютс€ NaN
end

MAX = cell(n_max,1);
for j = 1:n_max
    for i = 2 : 1 : size(temp{j},2)-1
        [M,I] = max([ temp{j}(i-1), temp{j}(i), temp{j}(i+1) ]);
        if I == 2
            MAX{j} = cat(1, MAX{j}, [find( Px(j,:) == M ), M]);
        end
    end
end

for i = 1:n_max
  if isempty(MAX{i})
      MAX{i} = [1, NaN];
  end
end

end

