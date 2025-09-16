function [im ,newDIM] = icatb_imageResize(image,DIM)
%  [im ,newDIM] = icatb_imageResize(image,DIM)
%  --------------------------------------------------------------------------
%   CREATED BY: Eric Egolf 
%   LAST MODIFIED: 12-29-03
%   ABOUT: Interpolates a 3d image to the size specified
%  
%   #########################################################################
%   
%   USAGE:
%  
%   - [im ,newDIM] = icatb_imageResize(image,DIM)
%     INFO: Interpolates a 3d image to the size specified
%     PARAMETER:
%        image = 3d(x,y,z) image
%        DIM = desired dimensions of the image
%
%     OUTPUT:
%        im = resized image
%        newDIM = image's new dimensions
%   #############################################################
%   
%   LICENSING:
%   
%  ------------------------------------------------------
%CODE FOR CHANGING Z DIMENSION

for i = 1:size(image,1)
    if(DIM(3)>1)
        im(i,:,:) = icatb_imgrescale( reshape(image(i,:,:),size(image,2),size(image,3)) ,[size(image,2),DIM(3)]);        
    else
        im(i,:,:) = reshape(image(i,:,:),size(image,2),size(image,3));    
    end
end
   

%CODE FOR CHANGING X and Y DIMENSION
if(DIM(2)>1 & DIM(1)>1)
    for i=1:size(im,3)
        im2(:,:,i) = icatb_imgrescale( reshape(im(:,:,i),size(im,1),size(im,2)) , [DIM(1),DIM(2)] );
    end
end

im=im2;
newDIM(1)=size(im,1);
newDIM(2)=size(im,2);
newDIM(3)=size(im,3);
