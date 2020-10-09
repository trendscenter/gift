function cm = icatb_getColormap(numOfComp,imageValues,useStructural, displayType)
% cm = icatb_getColormap(numOfComp,negValues,useStructural)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: gets colormap based on number of components, whether the data is
%   absoulte value and if there is a structural image. A matlab file specified 
%   by icatb_defaults contains the colorsmaps it uses
%
% #########################################################################
% 
% USAGE:
%
% -cm = icatb_getColormap(numOfComp,negValues,useStructural)
%   INFO: Gets colormap
%   ARGUMENTS:
%       numOfComp = number of components you need a color map for
%       imageValues = 1 pos neg
%                       2 pos
%                       3 abs
%       useStructural = 1 if you are overlaying on a structural image
%
%   OUTPUT:
%       cm = colormap
%
% #############################################################
% 
% LICENSING:
% 
% 
%------------------------------------------------------




%load colorfile defaults
icatb_defaults;
global COLORMAP_FILE;

if ~exist('displayType', 'var')
    displayType = 'other';
end

load(COLORMAP_FILE);

%--set up colormap
colorLength = 64;
for i=1:numOfComp
    %get colormaps
    if(i ==1)        
        if strcmp(lower(displayType), 'composite')
            cm = redV;
        else            
            if(imageValues == 2 | imageValues ==3)
                cm = hotN;
            elseif(imageValues == 4) % colormap for negative values
                cm = cold;
                % flip the color range for negative image values
                [nrows, ncols] = size(cold);                
                for ii = 1: nrows
                    cm(ii, 1:ncols) = cold(nrows - ii + 1, 1:ncols);
                end
            else
                cm = coldhot2; %coldhot;
            end
        end
    elseif(i==2)
        if(imageValues == 2 | imageValues ==3)
            cm = blueV; %green;
        else
            cm = blueV; %green;
        end
    elseif(i==3)
        if(imageValues == 2 | imageValues ==3)
            cm = greenV; %red;
        else
            cm = greenV; %red;
        end
    elseif(i==4)
        if(imageValues == 2 | imageValues ==3)
            cm = pinkV; %purple;
        else
            cm = pinkV; %purple;
        end
    elseif(i==5)
        if(imageValues == 2 | imageValues ==3)
            cm = yellowV; %blue;
        else
            cm = yellowV; %blue;
        end
    else
        if(imageValues == 2 | imageValues ==3)
            cm = yellowV; %orange; 
        else
            cm = yellowV; %orange; 
        end
    end
    
    colorbarSkip = size(cm,1)/colorLength;
    if(i==1)
        tempCM = [cm(1:colorbarSkip:end,:)];
    else
        tempCM = [tempCM;cm(1:colorbarSkip:end,:)];
    end
end
cm = tempCM;

if (useStructural)
    gm = linspace(0, 1, 64);
    if size(gm, 1) == 1
        gm = gm';
    end
    structCM = [gm, gm, gm];
    %structCM = gray;
    colorbarSkip = size(structCM,1)/colorLength;
    cm=[cm;structCM(1:colorbarSkip:end, :)];
end

colorbarSkip = ceil(size(cm,1)/256);
cm = cm(1:colorbarSkip:end, :);