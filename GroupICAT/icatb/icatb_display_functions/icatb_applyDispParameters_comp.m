function [icasig, HInfo] = icatb_applyDispParameters_comp(icasig, convertToZ, returnValue, threshValue, ...
    structDIM, HInfo)
% apply display parameters to the component images


% convertToZScores
if convertToZ
    icasig = icatb_convertImageToZScores(icasig);
end

%get desired image values
if returnValue == 1
    icasig = icasig;
elseif returnValue == 2
    temp = zeros(size(icasig));
    indices = find(icasig>0);
    temp(indices) = icasig(indices);
    icasig = temp;
    clear temp;
elseif returnValue == 3
    icasig = abs(icasig);

    % negative image values
elseif returnValue == 4
    temp = zeros(size(icasig));
    indices = find(icasig < 0);
    temp(indices) = icasig(indices);
    icasig = temp;
    clear temp;
end

%threshold image
for i=1:size(icasig,1)
    temp = zeros(size(icasig(i,:)));
    indices= find(abs(icasig(i,:)) >= threshValue);
    temp(indices) = icasig(i,indices);
    icasig(i,:) = temp;
end
