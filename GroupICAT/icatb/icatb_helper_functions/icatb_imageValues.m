function [promptString, matchIndex] = icatb_imageValues(index)

promptString = str2mat('Positive and Negative', 'Positive', 'Absolute', 'Negative');

indices = (1:size(promptString, 1));

if nargin > 1

    error('Accepts only one input argument');
    
end

matchIndex = 1;

if nargin == 1

    if ischar(index)

        matchIndex = strmatch(lower(index), lower(promptString), 'exact');
        if isempty(matchIndex)
            error('Check the input argument passed to image values');
        end

    elseif isnumeric(index)

        if round(index) ~= index
            error('Index must be an integer');
        end

        if index < 0 | index > size(promptString, 1)
            error('Check the index parameter passed');
        end

        matchIndex = index;

    else

        error('Index must be an integer or character');

    end

end