function icatb_generate_tal(cell_array, outname, voxsize, volStr, valStr)

gyri = cell2struct(cell_array{1,6},'name',2);
brodmann = cell2struct(cell_array{1,8},'name',2);
hemisphere = cell2struct(cell_array{1,4},'name',2);
t_values = cell_array{1,9};
x_value = cell_array{1,1};
y_value = cell_array{1,2};
z_value = cell_array{1,3};
mni_vals = cell_array{1, 10};
mni_vals = strtrim(regexprep(lower(mni_vals), 'mni:', ''));

%get unique gyri - This is asking for the labels that refer to things like
%superior temporal gyrus, cingulate gyrus, etc...
ind = 1;
clear gyrus;
gyrus(ind).name = gyri(1).name;
gyrus(1).leftcount = 0;
gyrus(1).rightcount = 0;
gyrus(1).BA = '';

for j = 1:size(gyri,1),
    flag = 0;
    for k = 1:length(gyrus),
        if (strcmp(gyrus(k).name,gyri(j).name)||isempty(gyri(j).name)),
            flag = 1;%matched a previously found gyrus
        end;
    end;
    if (~flag),%didn't match a previously found gyrus
        ind = ind+1;
        gyrus(ind).name = gyri(j).name;
        gyrus(ind).leftcount = 0;gyrus(ind).rightcount = 0;
        gyrus(ind).BA = '';
    end;
end;


for j = 1:size(gyri,1),
    for k = 1:length(gyrus),
        %which gyrus?
        if ((j==1)&&(k==1)),
            gyrus(k).RmaxT = [];
            gyrus(k).LmaxT = [];
        end;
        if (strcmp(gyrus(k).name,gyri(j).name)),
            if (isempty(strfind(gyrus(k).BA,brodmann(j).name))),
                gyrus(k).BA = strcat(gyrus(k).BA,':',brodmann(j).name);
            end;
            
            %             %get L/R
            if (strcmp(hemisphere(j).name,'Left Cerebrum')||strcmp(hemisphere(j).name,'Left Cerebellum')),
                gyrus(k).leftcount = gyrus(k).leftcount+1;
                if (gyrus(k).leftcount == 1), %get position and max T
                    gyrus(k).LmaxT = t_values(j);
                    gyrus(k).LmaxXYZ = [x_value(j) y_value(j) z_value(j)];
                    mniCoords = strread(mni_vals{j}, '%f');
                    gyrus(k).LMNImaxXYZ = [mniCoords(1), mniCoords(2), mniCoords(3)];
                end;
            elseif (strcmp(hemisphere(j).name,'Right Cerebrum')||strcmp(hemisphere(j).name,'Right Cerebellum')),
                gyrus(k).rightcount = gyrus(k).rightcount+1;
                if (gyrus(k).rightcount == 1), %get position and max T
                    gyrus(k).RmaxT = t_values(j);
                    gyrus(k).RmaxXYZ = [x_value(j) y_value(j) z_value(j)];
                    mniCoords = strread(mni_vals{j}, '%f');
                    gyrus(k).RMNImaxXYZ = [mniCoords(1), mniCoords(2), mniCoords(3)];
                end;
            end;
            
        end;
    end;
end;
for j = 1:length(gyrus),
    if (isempty(gyrus(j).LmaxT)),
        gyrus(j).LmaxT = -999;
        gyrus(j).LmaxXYZ = [0 0 0];
        gyrus(j).LMNImaxXYZ = [0, 0, 0];
    end;
    if (isempty(gyrus(j).RmaxT)),
        gyrus(j).RmaxT = -999;
        gyrus(j).RmaxXYZ = [0 0 0];
        gyrus(j).RMNImaxXYZ = [0, 0, 0];
    end;
end;

%gyrus(1).name = 'Brain Regions of Interest';

ba = formatBrodmann(cellstr(str2mat(gyrus.BA)));

clear a;
fid = fopen(outname, 'w+');

fprintf(fid, ['Area\t Brodmann Area\t', volStr, '\t', valStr, '\tMNI (x, y, z)\n']);

for j = 1:length(gyrus)
    fprintf(fid, '%s\t%s\t%3.1f/%3.1f\t%3.1f (%d, %d, %d)/%3.1f (%d, %d, %d)\t(%d, %d, %d)/(%d, %d, %d)\n', gyrus(j).name, ba{j}, gyrus(j).leftcount*voxsize, ...
        gyrus(j).rightcount*voxsize, gyrus(j).LmaxT, gyrus(j).LmaxXYZ(1), gyrus(j).LmaxXYZ(2), gyrus(j).LmaxXYZ(3), gyrus(j).RmaxT, ...
        gyrus(j).RmaxXYZ(1), gyrus(j).RmaxXYZ(2), gyrus(j).RmaxXYZ(3),  gyrus(j).LMNImaxXYZ(1),  gyrus(j).LMNImaxXYZ(2), gyrus(j).LMNImaxXYZ(3), ...
        gyrus(j).RMNImaxXYZ(1),  gyrus(j).RMNImaxXYZ(2), gyrus(j).RMNImaxXYZ(3));
end

fclose(fid);


function out = formatBrodmann(ba)
% Format brodmann input string

% Match regular expression
[startInd, endInd, tokens] = regexp(ba, '\d+');

ind = find(icatb_good_cells(startInd) ~= 0);

out = repmat({'*'}, length(ba), 1);
if ~isempty(ind)
    % Loop over selected indices
    for nInd = ind
        
        % Starting indices
        startv = startInd{nInd};
        % Ending indices
        endv = endInd{nInd};
        d = zeros(length(startv), 1);
        % Loop over starting vector
        for nd = 1:length(startv)
            d(nd) = str2num(ba{nInd}(startv(nd):endv(nd)));
        end
        
        d = unique(d);
        
        str = [];
        for nd = 1:length(d)
            tempStr = num2str(d(nd));
            if (nd == 1)
                str = tempStr;
            else
                str = [str, ', ', tempStr];
            end
        end
        
        % End loop over starting vector
        out{nInd} = str;
    end
    % End loop over selected indices
end




