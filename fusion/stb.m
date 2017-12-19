function out = stb(in, nMLAs, nOverlappingMLAs)

if nMLAs == 1
    out = in;
    return;
end

keepEdgeBeams = 1;

[nSamples, nBeamsIn, nPlanes] = size(in);

nTx           = nBeamsIn / nMLAs;
nBeamsEdge    = 2 * nOverlappingMLAs;
nBeamsKept    = 2 * nTx * (nMLAs/2 - nOverlappingMLAs);
nBeamsOverlap = (nTx-1) * nOverlappingMLAs;
nBeamsOut     = nBeamsKept + nBeamsOverlap + nBeamsEdge;

% outData = inData.prepareNewDataSet('Udt', zeros(nSamples,nBeamsOut,nPlanes), inData.params);
out = zeros(nSamples,nBeamsOut,nPlanes);

% Edge beams
edge_ix      = [1:(nBeamsEdge/2), (nBeamsIn+1-nBeamsEdge/2):nBeamsIn];
edge_out_ix  = [1:(nBeamsEdge/2), (nBeamsOut+1-nBeamsEdge/2):nBeamsOut];
out(:,edge_out_ix,:) = in(:,edge_ix,:);

% Kept beams
kept_ix     = zeros(1,nBeamsKept);
kept_out_ix = zeros(1,nBeamsKept);
nKeptPerTx  = nMLAs - 2 * nOverlappingMLAs;
for tx = 1:nTx
    kept_ix(     (tx-1)*nKeptPerTx + (1:nKeptPerTx) ) = (tx-1)*nMLAs + nOverlappingMLAs + (1:nKeptPerTx);
    kept_out_ix( (tx-1)*nKeptPerTx + (1:nKeptPerTx) ) = (tx-1)*(nMLAs-nOverlappingMLAs) + nOverlappingMLAs + (1:nKeptPerTx);
    
end
out(:,kept_out_ix(:),:) = in(:,kept_ix(:),:);

% Merged beams
stb_ix      = zeros(2,nBeamsOverlap);
stb_out_ix  = zeros(1,nBeamsOverlap);
for tx = 2:nTx
    stb_ix( 1,  (tx-2)*nOverlappingMLAs + (1:nOverlappingMLAs) ) = (tx-1)*nMLAs - nOverlappingMLAs + (1:nOverlappingMLAs);
    stb_ix( 2,  (tx-2)*nOverlappingMLAs + (1:nOverlappingMLAs) ) = (tx-1)*nMLAs + (1:nOverlappingMLAs);
    stb_out_ix( (tx-2)*nOverlappingMLAs + (1:nOverlappingMLAs) ) = (tx-1)*(nMLAs-nOverlappingMLAs) + (1:nOverlappingMLAs) - cSound.Utils.iif(keepEdgeBeams, 0, nOverlappingMLAs);
end
w  = linspace(1/(2*nOverlappingMLAs),1,2*nOverlappingMLAs);
w1 = fliplr(w(1:2:end));
w2 = w(1:2:end);
w1allBeams = repmat(w1, 1, nBeamsOverlap/nOverlappingMLAs);
w2allBeams = repmat(w2, 1, nBeamsOverlap/nOverlappingMLAs);

w1allBeams = repmat(w1allBeams(1,1:nBeamsOverlap), [nSamples, 1, nPlanes]);
w2allBeams = repmat(w2allBeams(1,1:nBeamsOverlap), [nSamples, 1, nPlanes]);

R = 0;
out(:,stb_out_ix(:),:) = in(:,stb_ix(1,:),:) .* w1allBeams .* exp(1j.*w2allBeams.*angle(R)) ...
    + in(:,stb_ix(2,:),:) .* conj(w2allBeams .* exp(1j.*w1allBeams.*angle(R)));

end