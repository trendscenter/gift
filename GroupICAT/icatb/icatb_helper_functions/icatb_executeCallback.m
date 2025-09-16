function icatb_executeCallback(objHandle)
% Execute the function callback using feval function as it takes into
% account the sub functions and any function overloading mechanism

% get the callback property for the object
objCallback = get(objHandle, 'callback');

if ~isempty(objCallback)
    
    % Execute the callback using feval
    if length(objCallback) > 1
        feval(objCallback{1}, objHandle, [], objCallback{2:end});    
    else    
        feval(objCallback{1}, objHandle);    
    end
    
end