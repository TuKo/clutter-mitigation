function writeMeasurementWindow(dataset, InX, InY, OutX, OutY)

if ~exist(dataset,'file')
    [success, ~,~ ] = mkdir(dataset);
    if ~success
        error('Cannot create folder');
    end
end

%read info file
inifile([dataset '/window.txt'],'write', ...
                  {'window','','In', [min(InX) max(InX) min(InY) max(InY)]; 'window','','Out', [min(OutX) max(OutX) min(OutY) max(OutY)]});

fprintf('In size: %d by %d pixels, area %d pixels\n', abs(InX(1)-InX(2))+1, abs(InY(1)-InY(2))+1, (abs(InY(1)-InY(2))+1)*(abs(InX(1)-InX(2))+1));
fprintf('In size: %d by %d pixels, area %d pixels\n', abs(OutX(1)-OutX(2))+1, abs(OutY(1)-OutY(2))+1, (abs(OutY(1)-OutY(2))+1)*(abs(OutX(1)-OutX(2))+1));
end
