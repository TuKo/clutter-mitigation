function Dout = getPhase(D, N, M)

dN = ceil(sqrt(size(D,2)));
% dM = ceil(size(D,2)/dN);
Dout = zeros(size(D));
for i = 1:size(D,2);
    atom = reshape(angle(D(:,i)),[N M]);
    atom = unwrap(atom,[],1);   
    Dout(:,i) = atom(:);
end

end