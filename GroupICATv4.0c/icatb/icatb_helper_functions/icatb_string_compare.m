function result = icatb_string_compare(strEntered, match)
%% Return true or false if the pattern exists in the string
%

strEntered = lower(strEntered);
match = lower (match);

chk = icatb_findstr(strEntered, match);

result = true;
if (isempty(chk))
    result = false;
end