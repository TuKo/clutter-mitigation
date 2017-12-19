function atoms = findClutterAtoms(D, THR, method, M, N, P)
% FINDCLUTTERATOMS returns a boolean vector describing for each atom in 
% dictionary D if it represents clutter or not.
%
% FINDCLUTTERATOMS(D, THR, METHOD)
% D - A row normalized dictionary.
% THR - The threshold value to decide whether an atom represents clutter
% METHOD - Options:
%                   1 - SVF/SVD (Using the first singular value)
%                   2 - TV measure in time (horizontal only).
% M, N, P - Patch/Atom size (M: axial, N: time, P: radial).

% Javier Turek
% GE Healthcare
% 26/8/2012

m = size(D,2);

valueAtom = zeros(m,1);
for k = 1:m
    atom = reshape(D(:,k),[M*P, N]);
    if method == 1
        [~,s,~] = svd(atom,0);
        s = diag(s)./sum(diag(s));
        valueAtom(k) = s(1);
    elseif method == 2
        valueAtom(k) = TVvalue(atom);
    else
        error('Wrong method selected');
    end
end

if method == 1
    atoms = (valueAtom > THR);
elseif method == 2
    valueAtom = valueAtom ./ max(valueAtom);
    atoms = (valueAtom < THR);
end

end