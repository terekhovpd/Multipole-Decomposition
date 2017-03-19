%ver 5.1
%19.03.2017
% Internal function for scattering patterns

function [ output_matrix   ] = dot_product_array_vector( Array,n1,n2,n3 )


output_matrix(1,1)=Array(1,1,1)*n1+Array(1,1,2)*n2+Array(1,1,3)*n3;
output_matrix(1,2)=Array(1,2,1)*n1+Array(1,2,2)*n2+Array(1,2,3)*n3;
output_matrix(1,3)=Array(1,3,1)*n1+Array(1,3,2)*n2+Array(1,3,3)*n3;


output_matrix(2,1)=Array(2,1,1)*n1+Array(2,1,2)*n2+Array(2,1,3)*n3;
output_matrix(2,2)=Array(2,2,1)*n1+Array(2,2,2)*n2+Array(2,2,3)*n3;
output_matrix(2,3)=Array(2,3,1)*n1+Array(2,3,2)*n2+Array(2,3,3)*n3;


output_matrix(3,1)=Array(3,1,1)*n1+Array(3,1,2)*n2+Array(3,1,3)*n3;
output_matrix(3,2)=Array(3,2,1)*n1+Array(3,2,2)*n2+Array(3,2,3)*n3;
output_matrix(3,3)=Array(3,3,1)*n1+Array(3,3,2)*n2+Array(3,3,3)*n3;







end
