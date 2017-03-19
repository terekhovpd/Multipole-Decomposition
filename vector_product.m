%ver 5.1
%19.03.2017
% Internal function for scattering patterns

function [ output_vec  ] = vector_product( a1,a2,a3,b1,b2,b3 )

output_vec(1)=    ( a2*b3 - a3*b2 );

output_vec(2)=    ( a3*b1 - a1*b3);

output_vec(3)=    ( a1*b2 - a2*b1 );

end

