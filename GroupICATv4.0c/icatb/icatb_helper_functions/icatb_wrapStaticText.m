function [staticTextH] = icatb_wrapStaticText(staticTextH)
% Input is static text handle
% output is the same but contains the wrapped text with the new string

inputStr = get(staticTextH, 'string');

% convert to cell
if ~iscell(inputStr)
    inputStr = {inputStr};
end

% wrap the text
[outputTextStr, newPos_static_text] = textwrap(staticTextH, inputStr);
% get the old position;
oldPos_staticText = get(staticTextH, 'position');
% set the new height
oldPos_staticText(4) = newPos_static_text(4);
% set the wrapped string
set(staticTextH, 'string', outputTextStr, 'position', oldPos_staticText);
