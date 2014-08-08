% *****************************************************************************
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
% 
% Copyright (c) 2014 Marco Aurelio Barbosa Fagnani Gomes Lotz (marcolotz.com)
%
% The source code in this document is licensed under Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. You must 
% credit the author of the source code in the way specified by the author or
% licenser (but not in a way to suggest that the author or licenser has given 
% you allowance to you or to your use of the source code). If you modify,
% transform or create using this source code as basis, you can only distribute
% the new source code under the same license or a similar license to this one.

% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.

% To see a copy of the license, access:
% creativecommons.org/licenses/by-nc-sa/4.0/legalcode
% *****************************************************************************




% Lung Nodules Candidates Detector
% Author: Marco Aurélio Lotz
% v. 0.2 23/06/2014

clear all;
close all;

%% 1. Pre-Processing
    % I guess it's better to work with raw dicom image
    % without any kind of pre-processing.
    
    % Load the DICOM Image:
    [image] = dicomread('I1400001');

    figure (1)
    imshow(image,[]);
    title('Original DICOM image');
    
%% 2. Lung Extraction
    % According to the paper, set the gray threshold to 375. All values
    % below are considered either background or lung parenchyma.
    
    thresholdHU = -375;
    % Reference for this convertion:
    % http://imagej.1557.x6.nabble.com/Hounsfield-Unit-conversion-td3686965.html
    
    thresholdGray = thresholdHU + 1000;
    
    % Generates the first mask using Gray Threshold.
    mask = zeros(size(image,1),size(image,2));
    mask(find(image < thresholdGray)) = 1;
    
    % Change from double
    mask = uint16(mask);
    
    figure (2)
    imshow(mask,[]);
    title('Binary mask');
    
    % It is recommended to use a scanning ball (2D disk) of 50X50 pixels by
    % the CT paper. In this case I am using 5, since 50x50 looks too large.
    
    % Note, a larger scanning ball (30x30) allows the detections of juxta
    % pleural nodules. Although it merges the two lungs.
    % Check this feature on image I7200000
    sDiskSize = 5;
    se = strel('disk',sDiskSize);
    closedMask = imclose(mask,se);
    
    figure (3)
    imshow(closedMask,[]);
    title('Closed Binary Mask');
       
    % Now one needs to remove the background. In order to perform this, one
    % has to find the connected components of the image and remove the ones
    % that touch the margin.
    
    connectedComp = bwconncomp(closedMask);
    
    ccMask = closedMask;
    
    % Removes the components that touch the borders
    ccMask = imclearborder(ccMask);
    
    figure(4)
    imshow(ccMask,[]);
    title('Mask of Connected Components');
    
    % Fills any holes inside the lungs mask
    noHolesMask = imfill(ccMask,'holes');
    figure (5)
    imshow(noHolesMask,[]);
    title('Filled Holes');
    
    % Extracted Background and Lung Parenchyma
    extracted1 = noHolesMask.*image;
    figure (6)
    imshow(extracted1,[]);
    title('Threshold Extraction');
    
    % Removing Aorta. As one may see, the aorta has a really small
    % frequency when compared to the remaining of the lung, so the otsu
    % threshold should be enought to remove it. The others regions that may
    % be removed from the image may be restaured using opennings.

    % otsuLevel = graythresh(extracted1);
    % otsuMask = im2bw(extracted1,otsuLevel);
    
    % figure (6)
    % imshow(otsuMask,[]);
    % title('otsuMask');
    
    
    % Final image of this section:
    output = extracted1;

%% 3. Nodes Candidates Detection
    
    input = extracted1;
    
    % According to the paper:
    % 2003_LUNG_Pulmonary Nodule Detection using chest CT images
    % "A gray value between 110 and 120 is used to extract the nodule
    % candidates bigger than 5x5 pixels." The paper, however, uses a byte
    % per pixel, when the standard DICOM format uses 2 bytes per pixel
    % (uint16). Thus, the value needs to be converted.
    % TODO: Check the info above.
    
    % According to the paper:
    minimumRangeByte = 110;
    maximumRangeByte = 120;
    byteRange = 255;
    
    % Found this value using a rule of three, one needs to double check
    % how valid is this convertion
    maximumRangeDICOM = 2235;
    minimumRangeDICOM = minimumRangeByte * maximumRangeDICOM /byteRange;
    maximumRangeDICOM = maximumRangeByte * maximumRangeDICOM /byteRange;
    
    % Let's find the regions that are within this grey level:
    maskRegion = zeros(size(input,1),size(input,2));
    maskRegion(find(input > minimumRangeDICOM & input < maximumRangeDICOM)) = 1;
    
    figure (7)
    imshow(maskRegion,[]);
    title('Mask for regions in the grey boundary');
    
    % Selects only connected regions larger than 5x5 pixels:
    % (Connectivity 8)
    connectedCandidates = bwconncomp(maskRegion);
    
    % Bouding Box will be used to calculate the size of the component.
    CandidatesProps = regionprops(connectedCandidates, 'BoundingBox');
    
    numberOfCandidates = connectedCandidates.NumObjects;
    
    % New mask that will only have objects larger than 5x5 pixels.
    largeCandidatesMask = maskRegion;
    
    for count = 1:numberOfCandidates
        Buffer = CandidatesProps(count).BoundingBox;
        % If the bounding box is smaller than 5x5
        if (Buffer(3)<5 || Buffer(4) < 5)
            %Erases the candidates that fail the test.
            largeCandidatesMask(connectedCandidates.PixelIdxList{count}) = 0;
        end
    end
    
    figure (8)
    imshow(largeCandidatesMask,[]);
    title('Mask for candidates larger than 5x5 pixels');
