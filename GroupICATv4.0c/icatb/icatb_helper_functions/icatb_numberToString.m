function str = icatb_numberToString(nums)
% number to string. pad space after the number

try
    str = arrayfun(@num2str, nums, 'uniformoutput', false);
catch
    str = cell(length(nums), 1);
    for n = 1:length(nums)
        str{n} = num2str(nums(n));
    end
end

str = char(str);