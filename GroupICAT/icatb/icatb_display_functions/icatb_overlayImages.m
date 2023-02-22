function [im,maxICAIM,minICAIM,minInterval,maxInterval] = icatb_overlayImages(icasig,structuralImage,icaDIM,structDIM,imageValues)
%  [im,maxICAIM,minICAIM,minInterval,maxInterval] = icatb_overlayImages(icasig,structuralImage,icaDIM,structDIM)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: Puts components on top of structural image 
%
% #########################################################################
% 
% USAGE:
%  
%  -[im,maxICAIM,minICAIM,minInterval,maxInterval] = icatb_overlayImages(icasig,structuralImage,icaDIM,structDIM)
%   INFO: Puts components on top of structural image and scales the images
%   as needed
%   PARAMETERS:
%       icasig = 4d(number of components,x,y,z) spatial components
%       structuralImage = 3d (x,y,z) structural image
%       icaDIM = = dimensions of spatial components
%       structDIM = dimensions of structural components
%
%   OUTPUT:
%       im = image with components overlaid on structural
%       maxICAIM = vector of max ica values
%               ex: values the maxInterval for a particular component
%               represents. So if  maxInterval(component1) = 100 and my
%               maxICAIM = 5 then the value 100 in im really represents 5
%       minICAIM = same as maxICAIM but min
%       maxInterval = vector of maxInterval values
%               intervals that each component and the structural image is
%               scaled between. 
%               *if you had 3 components and 1 structural image component 1
%               would be scaled between 0 and 100 component two would be
%               scaled between 100+unitColor and 200, component 3 would be
%               scaled between 200+unitColor and 300, and the structural
%               image would be scaled between 300+unitColor and 400.
%               To find would what the scaling values really stand for use
%               maxICAIM and minICAIM.
%       minInterval = same as maxInterval but min
%       minValue = ?
%       maxValue = ?
% #############################################################
% 
% LICENSING:
% 
% 
%------------------------------------------------------

if(~exist('imageValues','var'))
    imageValues = 1;
end

%icasig = 4D
%struct = 4D
minInterval = 0;
maxInterval = 100;

%reshape structural
structFlat = reshape(structuralImage,structDIM(1)*structDIM(2)*structDIM(3),1);
clear structuralImage;
icaFlat=reshape(icasig,size(icasig,1),icaDIM(1)*icaDIM(2)*icaDIM(3));

numToOverlay = size(icasig,1);


for i=1:numToOverlay
        [tempMax] = max(icasig(i,:));
        [tempMin] = min(icasig(i,:));
        
        %case: When Empty Icasig
        if(tempMax==0 & tempMin ==0)
            icasig(i,1) = 1;
        end
        
        %case: Positive and Negative    
        if(imageValues ==1)
            
            if(tempMax>abs(tempMin))
                maxICAIM(i) = tempMax;
                minICAIM(i) = -tempMax;    
            else
                maxICAIM(i) = abs(tempMin);
                minICAIM(i) = tempMin;    
            end
        
        %case: Positive Only
        elseif(imageValues ==2)
            
            maxICAIM(i) = tempMax;
            minICAIM(i) = tempMin;
        
        %case: Absolute Value
        elseif(imageValues==3)
            maxICAIM(i) = tempMax;
            minICAIM(i) = tempMin;
            
        % case: Negative Value
        elseif(imageValues==4)
            maxICAIM(i) = tempMax;
            minICAIM(i) = tempMin; 
            
        %case: Error
        else
            %infoCell{1} = imageValues;
            %icatb_error('Invalid Image Value Parameter',InfoCell);    
            error('Invalid image value parameter');
        end
     
end
%   %---------------------Get maxs and mins for scaling
%   for i =1:numToOverlay  
%     minICAVal(i) = min(icaFlat(i,:));
%     maxICAVal(i) = max(icaFlat(i,:));
%     %--case when nothing is in image
%     if(minICAVal(i)==0 & maxICAVal(i)==0)
%         icaFlat(i,1) = 1;
%         maxICAVal(i)=1;
%     end    
%     minICAIM(i) = minICAVal(i);
%     maxICAIM(i) = maxICAVal(i);
%     
%     if(min(min(min(icasig)))<0)
%         negValues = 1;
%     else
%         negValues =0;
%     end   
%     if( abs(minICAVal(i)) > abs(maxICAVal(i)) & negValues)
%         maxICAIM(i) = -minICAVal(i);
%     elseif(negValues)
%         minICAIM(i) = -maxICAVal(i);
%     end
% end


%--Loop through each component to scale and overlay
for i =1:numToOverlay  
    %get unit color
    unitColor = icatb_range(icaFlat(i,:))/(min([256 64*(numToOverlay+1)])/(numToOverlay+1));
    
    indices(i).ind = find(icaFlat(i,:) ~= 0);
    
    %% If maxICAIM(i) is equal to zero
    if(maxICAIM(i) == 0) 
        maxICAIM(i) = 10^-8;
    end
    
    icaFlat(i,find(icaFlat(i,:)==maxICAIM(i))) = icaFlat(i,find(icaFlat(i,:)==maxICAIM(i)))-unitColor;
    
    %scale values in component from minInterval to maxInterval
    multi = maxICAIM(i)/icatb_range([minICAIM(i) maxICAIM(i)]);
    icaFlat(i,:) = ( (icaFlat(i,:)./maxICAIM(i)) +  abs(minICAIM(i))/maxICAIM(i)) * (multi) * maxInterval;
    icaFlat(i,:) = icaFlat(i,:) + (i-1)*maxInterval;
    
end


%scale structural
unitColor = icatb_range([minInterval maxInterval])*(numToOverlay+1)/min([256 64*(numToOverlay+1)]);
newMin = maxInterval*numToOverlay;
newMax = maxInterval*(1+numToOverlay);
structFlat = icatb_scaleImage(newMin,newMax,structFlat);



%overlay structural
overlayImage = structFlat(:);

%overlay components
for i =1:numToOverlay
    overlayImage(indices(i).ind) = icaFlat(i,indices(i).ind);
end
im=overlayImage;
overlayMap = [];