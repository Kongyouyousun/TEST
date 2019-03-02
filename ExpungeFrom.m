function C= ExpungeFrom( A, B)
%#codegen
 C=A(~ismember(A,B));
end
