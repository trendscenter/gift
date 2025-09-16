function [V, HInfo] = icatb_returnHInfo(P)
% Return Volume and Header information of a vector of 3D analyze images,
% one 4D nifti or one 4D analyze image
%
% Input:
%
% 1. P - file names as a character array (Pass either a nifti, 4D analyze
% file or a vector of 3D analyze files)
%
% Output:
%
% 1. V - volume information
% 2. HInfo - header information


if size(P, 1) > 1
    P = deblank(P(1, :));
end

gzipFile = 0;
if (icatb_findstr(lower(P), '.gz'))
    gzipFile = 1;
end

if (~gzipFile)
    if icatb_findstr(lower(P), '.gii')
        gd = gifti(P);
        V(1).fname = P;
        V(1).dim = [length(gd.cdata), 1, 1];
        VOX = [1, 1, 1];
        gd.cdata = [];
        V(1).gifti = gd;
    else
        V = icatb_get_vol_nifti(P);      
        if (V(1).private.hdr.dim(1) == 6)
            % treat it as cifti data
            V(1).cifti.dim = [V(1).private.hdr.dim(7), V(1).private.hdr.dim(6)];
            V(1).dim = [V(1).cifti.dim(1), 1, 1, V(1).cifti.dim(2)];
            dat = ft_read_cifti(icatb_parseExtn(P));
            dat.dtseries = [];
            V(1).cifti.dat=dat;
        end
        
        % voxel size
        VOX = double(V(1).private.hdr.pixdim(2:4));
    end
    
else
    [hdr, V] = icatb_read_gzip_nii(icatb_parseExtn(P), 'read_hdr_only', 1);
    VOX = hdr.pixdim(2:4);
end

% make the voxel sizes positive
VOX = abs(VOX);

%save Header info in structure
HInfo = struct('DIM', V(1).dim(1:3), 'V', V, 'VOX', VOX);