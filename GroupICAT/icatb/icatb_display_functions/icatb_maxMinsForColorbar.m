function [icasig,minValues,maxValues] =icatb_maxMinsForColorbar(imageValues,colorbarScale,icasig)
%[min,max] = function icatb_maxMinsForColorbar(imageValues,colorbarScale,icasig)
% imageValues = 1 pos and neg
%               2 = pos
%               3 = absolute
%
% icasig = (numberOfComp,x,y,z);

numOfComp = size(icasig,1);
xdim = size(icasig,2);
ydim= size(icasig,3);
zdim = size(icasig,4);
icasig=reshape(icasig,numOfComp,xdim*ydim*zdim);

for i=1:numOfComp
        [tempMax] = max(icasig(i,:));
        [tempMin] = min(icasig(i,:));
        
        %case: When Empty Icasig
        if(tempMax==0 & tempMin ==0)
            icasig(i,1) = 1;
            
        %case: Positive and Negative    
        elseif(imageValues ==1)
            
            if(tempMax>abs(tempMin))
                maxValues(i) = tempMax;
                minValues(i) = -tempMax;    
            else
                maxValues(i) = abs(tempMin);
                minValues(i) = tempMin;    
            end
        
        %case: Positive Only
        elseif(imageValues ==2)
            
            maxValues(i) = tempMax;
            minValues(i) = tempMin;
        
        %case: Absolute Value
        elseif(imageValues==3)
            maxValues(i) = tempMax;
            minValues(i) = tempMin;
            
        elseif(imageValues==4)
            maxValues(i) = 0;
            minValues(i) = tempMin;
            
        %case: Error
        else
            infoCell{1} = imageValues;
            icatb_error('Invalid Image Value Parameter',InfoCell);    
        end
    
    
    
end


% To Make Uniform
if(strcmp(colorbarScale ,'uniform'))
    maxValues(1:numOfComp) = max(maxValues);
    minValues(1:numOfComp) = min(minValues);
end


icasig=reshape(icasig,numOfComp,xdim,ydim,zdim);
