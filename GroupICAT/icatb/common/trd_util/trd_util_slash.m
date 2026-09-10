function [outputArg1,outputArg2] = trd_util_slash(s_path)
    % Creates a slash at the end of a path OS dependently
    
    if s_path(end) ~= filesep
        s_path = [s_path filesep];
    end        

end

