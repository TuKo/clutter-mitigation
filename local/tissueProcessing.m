function mov = tissueProcessing(movIQ, DR, GAIN)

% Tissue processing
% movIQ = correctMLA(movIQ);

% c =0.00025;
% mov = log(1+c*abs(movIQ));

if nargin < 2
    DR = 30;
    GAIN = -50;
elseif nargin < 3
    GAIN = -50;
end

mov = (20/DR)*(log10(1+abs(movIQ))+ GAIN/20);

end