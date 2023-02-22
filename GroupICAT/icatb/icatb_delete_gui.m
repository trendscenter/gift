function icatb_delete_gui(figTag)
% Delete guis with the tag specified

showHiddenHandles = get(0, 'ShowHiddenHandles');

try
    set(0, 'ShowHiddenHandles', 'on');
    
    childrens = get(0, 'Children');
    
    if ~iscell(figTag)
        figTag = {figTag};
    end
    
    % Loop over tags
    for nTag = 1:length(figTag)       
        ind = findobj(0, 'tag', figTag{nTag});
        if ~isempty(ind)
            delete(ind);
        end
    end
    % End loop over tags
    
    set(0, 'ShowHiddenHandles', showHiddenHandles);
    
catch
    set(0, 'ShowHiddenHandles', showHiddenHandles);  
end