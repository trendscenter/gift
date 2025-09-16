function [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, optional)
% name complex images

icatb_defaults;

% variables to detect the reading and writing of complex images
global READ_NAMING_COMPLEX_IMAGES;
global WRITE_NAMING_COMPLEX_IMAGES;


complexInfo = [];

% naming of complex images
if ~isfield(sesInfo.userInput, 'read_complex_images')
    sesInfo.userInput.read_complex_images = 'real&imaginary';
end

if ~isfield(sesInfo.userInput, 'write_complex_images')
    sesInfo.userInput.write_complex_images = 'real&imaginary';
end

read_complex_images = sesInfo.userInput.read_complex_images;
write_complex_images = sesInfo.userInput.write_complex_images;

% file naming for reading complex images
if strcmpi(read_complex_images, 'real&imaginary')
    sesInfo.userInput.read_complex_file_naming = READ_NAMING_COMPLEX_IMAGES.real_imag;
else
    sesInfo.userInput.read_complex_file_naming = READ_NAMING_COMPLEX_IMAGES.mag_phase;
end

% file naming for writing complex images
if strcmpi(write_complex_images, 'real&imaginary')
    sesInfo.userInput.write_complex_file_naming = WRITE_NAMING_COMPLEX_IMAGES.real_imag;
else
    sesInfo.userInput.write_complex_file_naming = WRITE_NAMING_COMPLEX_IMAGES.mag_phase;
end

if strcmpi(optional, 'read')
    % complex file naming
    complexInfo.complex_file_naming = sesInfo.userInput.read_complex_file_naming;
    % complex type
    complexInfo.complexType =  sesInfo.userInput.read_complex_images;
else
    % complex file naming
    complexInfo.complex_file_naming = sesInfo.userInput.write_complex_file_naming;
    % complex type
    complexInfo.complexType =  sesInfo.userInput.write_complex_images;
end