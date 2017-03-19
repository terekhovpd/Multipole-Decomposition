%ver 5.1
%19.03.2017
% Internal function for scattering patterns

function [ output_vector  ] = dot_product_matrix_vector( M,n1,n2,n3 )


output_vector(1)=M(1,1)*n1+M(1,2)*n2+M(1,3)*n3;
output_vector(2)=M(2,1)*n1+M(2,2)*n2+M(2,3)*n3;
output_vector(3)=M(3,1)*n1+M(3,2)*n2+M(3,3)*n3;


end


