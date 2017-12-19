function [InX, InY, OutX, OutY] = readMeasurementWindow(dataset)

%read info file
[results]=inifile([dataset '/window.txt'],'read', ...
                  {'window','','In'; 'window','','Out'});

[In] = str2num(results{1});
[Out] = str2num(results{2});

[InY, InX] = meshgrid(In(3):In(4),In(1):In(2));
InX = InX(:);
InY = InY(:);

[OutY, OutX] = meshgrid(Out(3):Out(4),Out(1):Out(2));
OutX = OutX(:);
OutY = OutY(:);

end
