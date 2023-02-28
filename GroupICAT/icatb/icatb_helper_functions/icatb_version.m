function sVer = icatb_version()
% Fucntion that reads version number from README.md
% Function returns a string (sVer) which is the version of GIFT
% E.g., sVer = v4.0.4.0

sDirBefore = pwd; 

sVer = 'version not found'; % Just in case
[sPathHere sFile sExt] = fileparts(which('icatb_version.m')); %Filepath to this file
cd(sPathHere);
try
    fid = fopen('../../../README.md','rt');
    if fid == -1
        error('Cannot open %s.',filename);
    end
    M = struct('vertices',[],'faces',[]);
    %-Read File
    l = fgetl(fid); %Version is not supposed to be at first line
    while ~feof(fid)
        l = fgetl(fid);
        if strcmp(l,'<!-- PLEASE DO NOT EDIT THIS LINE OR LINE BELOW -->')
            l = fgetl(fid);
            sVer = l(37:end);
            break;
        end
    end
    cd(sDirBefore);
    s = fclose(fid);
catch
    cd(sDirBefore);
    s = fclose(fid);
end
