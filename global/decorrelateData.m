function IQdata = decorrelateData(data, rho, varargin)

if rho ~= 1.0
    R=[1 rho; rho 1];
    w=chol(R);
else
    w=[1 1; 0 0];
end
C1 = w(1,2);
C2 = w(2,2);

nFrames = size(data,3);

if nargin < 3
    ord = randintrlv(1:nFrames,sum(10000*clock));
else
    ord = varargin{1};
end
if nFrames < 2
    IQdata = data;
    return;
end

IQdata = zeros(size(data,1),size(data,2),size(data,3)-1);
IQdata(:,:,1) = C1*data(:,:,ord(1)) + C2*data(:,:,ord(2));

for k = 2:nFrames-1
    IQdata(:,:,k) = C1*IQdata(:,:,k-1) + C2*data(:,:,ord(k+1));
end

end