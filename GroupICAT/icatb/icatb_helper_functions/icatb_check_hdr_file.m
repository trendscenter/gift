function  icatb_check_hdr_file(file_to_check)      
% sub function to check if there is a header file or not for image file

hdrFile1 = [file_to_check(1:end-3), 'hdr'];
hdrFile2 = [file_to_check(1:end-3), 'HDR'];
if ~(exist(hdrFile1, 'file') || exist(hdrFile2, 'file'))
    error(['Header file for image ', file_to_check, ' doesn''t exist']);
end