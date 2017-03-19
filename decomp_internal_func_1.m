%Version 5.1
%19.03.2017
%Internal function for data_concatenation, allows to merge different frequency step for equal parameter value

function [ mask ] = decomp_internal_func_1 ( mask, i, number )


for j = 1:number
	mask (i,j) = max ( mask(i,j) , mask(i+1,j) );
end
mask (i+1,:) = [];

end
