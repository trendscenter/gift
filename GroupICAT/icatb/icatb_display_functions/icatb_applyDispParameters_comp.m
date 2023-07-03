function [icasig, HInfo] = icatb_applyDispParameters_comp(icasig, convertToZ, returnValue, threshValue, ...
    structDIM, HInfo)
% apply display parameters to the component images


% convertToZScores
if convertToZ
    icasig = icatb_convertImageToZScores(icasig);
end

icasig = applyDispParams(icasig, returnValue, threshValue);

% %get desired image values
% if returnValue == 1
%     icasig = icasig;
% elseif returnValue == 2
%     temp = zeros(size(icasig));
%     indices = find(icasig>0);
%     temp(indices) = icasig(indices);
%     icasig = temp;
%     clear temp;
% elseif returnValue == 3
%     icasig = abs(icasig);
% 
%     % negative image values
% elseif returnValue == 4
%     temp = zeros(size(icasig));
%     indices = find(icasig < 0);
%     temp(indices) = icasig(indices);
%     icasig = temp;
%     clear temp;
% end
% 
% %threshold image
% for i=1:size(icasig,1)
%     temp = zeros(size(icasig(i,:)));
%     indices= find(abs(icasig(i,:)) >= threshValue);
%     temp(indices) = icasig(i,indices);
%     icasig(i,:) = temp;
% end


function data = applyDispParams(data, returnValue, threshold)
%% Apply display parameters

threshold = abs(threshold);

if (returnValue == 2)
    % Positive
    data(data < 0) = 0;
elseif (returnValue == 3)
    % Absolute value
    data = abs(data);
elseif (returnValue == 4)
    % Negative value
    data(data > 0) = 0;
end

if (length(threshold) > 1)
    
    if ((returnValue == 2) || (returnValue == 3))
        % Positive or absolute
        data(data < min(threshold)) = 0;
        data(data > max(threshold)) = max(threshold);
    elseif (returnValue == 4)
        % Negative
        data(abs(data) < min(threshold)) = 0;
        data(abs(data) > max(threshold)) = -max(threshold);
    else
        % Positive and Negative
        neg_threshold = -threshold;
        
        pos_inds = (data > 0);
        neg_inds = (data < 0);
        
        % Handle positive range
        tmp1 = data(pos_inds);
        tmp1(tmp1 < min(threshold)) = 0;
        tmp1(tmp1 > max(threshold)) = max(threshold);
        
        % handle negative range
        tmp2 = data(neg_inds);
        tmp2(tmp2 > max(neg_threshold)) = 0;
        tmp2(tmp2 < min(neg_threshold)) = min(neg_threshold);
        
        data(pos_inds) = tmp1;
        data(neg_inds) = tmp2;
    end
    
    %data(abs(data) < min(threshold)) = 0;
    %data(abs(data) > max(threshold)) = 0;
else
    data(abs(data) < threshold) = 0;
end
