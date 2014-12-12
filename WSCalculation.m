%% WATERSHED SET OF FUNCTIONS ------
function partitioned_emg = WSCalculation(imi,RelativeEffortLevel)
%
% function partitioned_emg = emgpartition(imi);
%
%

imi(1) = (imi(1,2)+imi(2,1)+imi(2,2))/3;

% changing image contrast if specified
if nargin == 2 && ~isempty(RelativeEffortLevel) % then image contrast is enhanced
    imi = equalization(imi,RelativeEffortLevel);
end

% concatenating image frame border with nearest pixel values
A = [imi(1) imi(1,:) imi(1,end); imi(:,1) imi imi(:,end);imi(end,1) imi(end,:) imi(end)];

% creting sobel kernel
Kx = [1 2 1; 0 0 0;-1 -2 -1]; % horizontal edge detector
Ky = [1 0 -1; 2 0 -2; 1 0 -1]; % vertical edge detector

% evaluating partial derivatives along x (rows) and y (columns) directions
Gx = ifft2(fft2(A).*conj(fft2(Kx,size(A,1),size(A,2))));
Gy = ifft2(fft2(A).*conj(fft2(Ky,size(A,1),size(A,2))));

G = sqrt(Gx.^2 + Gy.^2); % the gradient of emg image
G = G(1:end-2,1:end-2); % removing the initially concatenated pixels

% smoothing the gradient with morphological operations "opening" and "closing"
% imamge opening removes noise on the extremity of objects in the image, as
% well as narrow connections, whereas image closing fills narrow long gulfs and holes 

% smooth_G = immorph(immorph(G,'erode'),'dilate'); % image opening
% smooth_G = immorph(immorph(smooth_G,'dilate'),'erode'); % image closing

smooth_G = immorph(immorph(G,'dilate'),'erode'); % image closing
smooth_G = immorph(immorph(smooth_G,'erode'),'dilate'); % image opening

% at last, the watershed algorithm is applied to segment the image
partitioned_emg = water(smooth_G);

for i = 1:max(partitioned_emg)
    m(i) = mean(abs(imi(partitioned_emg==i)));
end
[junk,clusters_with_max_amplitude] = max(m);
partitioned_emg(partitioned_emg~=clusters_with_max_amplitude) = 0;


function modified_im = immorph(imi,type,mask)
%
% function modified_im = immorph(imi)
%
% this funciton applies two basic operation of Mathematical Morphology to
% the scaled image "imi", which can be of any class.
% 
% image dilatation or erosion is applied according to the input argument
% "type":
% type = 'erode' | type = 'dilate' | type = 'denoising'
%
% denoising is a customized eroding operation to remove isoleted borders
% and narrow connections within the same group of channels, segmented with
% emgpartition function. In this case "imi" is an image with labelled group
% of pixels, and each group is individually scanned to remove the noise
%
% mask argument can be any structure element, defined as binary images
% of size [3 x 3]. (the algorithm may be improved to use any mask)
% 
%                                 0 0 0
% By default, the mask is set to  0 0 0 , which correspond to a flat
%                                 0 0 0
%
% structuring element appropriate for grayscale (scaled image)
% morphological operations
%
% This function is called by: "emgpartition"
%
% Author: Taian Vieira
% Date: 05/09/08

% concatenating image frame border with nearest pixel values
if strmatch(type,'denoising')
    A = zeros(size(imi)+2);
    A(2:end-1,2:end-1) = imi; % zeros are concatenated to the borders
    imi = A;
    
    for i = min(A(:)):max(A(:))
        if i==0, continue, end
        ind = find(A == i); % geting the index of pixels labelled as "i"
        
        for k = 1:numel(ind) % these pixels are scanned for ones(2,2) connectivity
            [x,y] = ind2sub(size(A),ind(k));
            neighbors = A([x-1 x x+1],[y-1 y y+1]);
            
            if ~( all(neighbors([1 2 4])==A(x,y)) | all(neighbors([2 3 6])==A(x,y))...
                    | all(neighbors([4 7 8])==A(x,y)) | all(neighbors([6 8 9])==A(x,y)) )
                imi(x,y) = 0;
            end
        end
    end
    modified_im = imi(2:end-1,2:end-1);
    return
end

if nargin == 2
    mask = zeros(3,3); % Flat mask for grayscale image
else
    mask = flipud(fliplr(mask));
end

mask_origin = (size(mask)+1)/2;
mask_offset = size(mask)-mask_origin;

A = [imi(1) imi(1,:) imi(1,end); imi(:,1) imi imi(:,end);imi(end,1) imi(end,:) imi(end)];
% [Nrows,Ncols] = size(imi);
% A = interp2(imi,linspace(1,Ncols,Ncols+2*mask_offset(2))',linspace(1,Nrows,Nrows+2*mask_offset(1)),'nearest');
% A([0:size(imi,1)-1]+mask_origin(1),[0:size(imi,2)-1]+mask_origin(2)) = imi;

modified_im = zeros(size(imi)); % initializing output image

% for row = [0:size(imi,1)-1]+mask_origin(1)
% 
%     for col = [0:size(imi,2)-1]+mask_origin(2)
%         I = mask + A([row-mask_offset(1)]:[row+mask_offset(1)],[col-mask_offset(2)]:[col+mask_offset(2)]);
%         
%         if strmatch(type,'erode')
%             modified_im([row-mask_offset(1)],[col-mask_offset(2)]) = min(I(:));
%     
%         elseif strmatch(type,'dilate')
%             modified_im([row-mask_offset(1)],[col-mask_offset(2)]) = max(I(:));
%             
%         end
%     end
% end
% 
for row = 2:size(A,1)-1

    for col = 2:size(A,2)-1
        I = mask + A(row-1:row+1,col-1:col+1);
        
        if strmatch(type,'erode')
            modified_im(row-1,col-1) = min(I(:));
    
        elseif strmatch(type,'dilate')
            modified_im(row-1,col-1) = max(I(:));
            
        end
    end
end


function seg_im = water(imi)
%
% function seg_im = water(imi)
%
% this function applies the watershed transform to the input image "imi".
% The watershed algorithm was proposed by (Vincent and Soille, 1991) and
% segments the image according to:
%
% Sorting  -> all pixels are sorted with increasing intensity
% 
% Flooding -> pixels are directly accessed due to the sorting procedure
% and compared for three different conditions. 1) If it is close to a
% single group of already labelled neighbors then it receives the label of
% the group. 2) If it has no neighbors then it is a new minimum and
% receives a new label; 3) If it has neighbors with different labels, then it
% is assigned the watershed label
%
% The resulting watershed line depends on the connectivity mask used for
% computing the neighbors. By default a 8-connectivity is used, thus
% producing watershed with 4-connected neighbors. 4-connectivity results in
% high selectivity for the correlation map of the sEMG signals (originally
% applied to assess load sharing between both gastrocnemius during standing)
%
% This function calls NGneighbors for returning the index of neighbor pixels
%
%  Reference: L. Vincent and P. Soille, "Watershed in digital spaces: 
%  An efficient algorithm based on immersion simulations," IEEE
%  Transactions on Pattern Analysis and Machine Intelligence, vol. 13, n. 6,
%  1990, pp. 583-598.
%
% Author: Taian Vieira
%   Date: 05/09/08

MASK = -2;
WSHED = 0;
FICTITIOUS = -1;
currentLabel = 0;

L = -ones(size(imi)); % initializing the whatershed output
    
% Initializing the pixel queue.
auxqueue = 1;

% Initializing the distance array with zeros
imd = zeros(size(imi));

[h,ind] = sort(imi(:)); % sorting the image with increasing order of pixels intensity
N = length(ind);

numProcessedPixels = 1;
while (numProcessedPixels <= N)

    % Find the next set of pixels with the same value.
    k1 = numProcessedPixels;
    currentLevel = h(k1);
    k2 = k1 + 1;

    while ((k2 < N) & (h(k2) == currentLevel)), k2 = k2+1;

             end % increment the pointer to the last pixel with the same "currentlevel" intensity

    k2 = k2-1; % correct for the last + 1 from the while loop

    % Mask all image pixels whose value equals currentLevel.
    for k = k1:k2
        p = ind(k); % index pointing to the pixel with intensity equal to "currentlevel"
        L(p) = MASK; % masking this pixel

        s=NGneighbors(imi,p); % geting the index of neighbor pixels 8-connectivity

        for i=1:length(s);
            q=s(i);
            if ((L(q) > 0) | (L(q) == WSHED)) % If it is labelled to a catchment basin or it is a watershed

                % Initialize queue with neighbors at currentLevel of current basins or watersheds.
                imd(p) = 1; % place 1 at the p position of the work image (distance image)
                queue(auxqueue) = p;
                auxqueue = auxqueue+1; % auxqueue corresponds to the ptr_last of soille algorithm

                break;
            end
        end
        numProcessedPixels = numProcessedPixels+1;
    end

    currentDistance = 1;
    queue(auxqueue) = FICTITIOUS; % fictitious index
    auxqueue = auxqueue+1; % incrementing ptr_last

% Extending the catchment basins.
    while true
        p = queue(1); % geting the first pixel stored in the queue
        queue(1) = []; % removing the first pixel from the queue
        auxqueue = auxqueue-1; % the ptr_last must be decremented since queue structure reduced

        if (p == FICTITIOUS)

            if isempty(queue) % checking if queue is empty
                break;

            else % there is still another sequence after the fictitious index
                queue(auxqueue) = FICTITIOUS; % fictitious index
                auxqueue = auxqueue+1; % incrementing ptr_last
                currentDistance = currentDistance+1;
                p = queue(1); % geting the first index of this new sequence
                queue(1) = [];
                auxqueue = auxqueue-1; % decrementing ptr_last
            end
        end

% search for already labelled neighbor pixels or watershed lines. The
% closes label is assigned to the new pixel if there is not watershed line
% in its neighborhood. Masked neighbors are stored into the queue and their
% distance is set to 1

        closestDist = 1000; % initializing closest distance with an arbitrarily higher number
        closestLabelValue = 0;
        IsUniqueclosestLabelValue = true;

        s=NGneighbors(imi,p);

        for i=1:length(s);
            q=s(i);

            if ((L(q) > 0) | (L(q) == WSHED))

                if (imd(q) < closestDist) % chiecking for the closest labelled pixel
                    closestDist = imd(q);

                    if (L(q) > 0)
                        closestLabelValue = L(q);
                    end

                elseif imd(q) == closestDist % there is another labeled pixel in the neighborhood

                    if (L(q) > 0)

                        if ((closestLabelValue > 0) & (L(q) ~= closestLabelValue)) % then, there are distinct groups on the neighborhood
                            IsUniqueclosestLabelValue = false;
                        end
                        closestLabelValue = L(q);
                    end
                end

            elseif ((L(q) == MASK) && (imd(q) == 0)) % q and p have same intensity and thus are plateau pixels
                imd(q) = currentDistance + 1; % the distance of p neighbors is longer from the q pixel
                queue(auxqueue) = q;
                auxqueue = auxqueue+1;
            end
        end

% At this step pixel p is labelled either as catchment basins or watershed line

        if ((closestDist < currentDistance) & (closestLabelValue > 0))

            if (IsUniqueclosestLabelValue & ((L(p) == MASK) | (L(p) == WSHED))) % then, p belongs to the closest basin
                L(p) = closestLabelValue;
                imd(p)=0;

            elseif ( ~IsUniqueclosestLabelValue | (L(p) ~= closestLabelValue)) % then, p is in between and thus is part of a wshed line
                L(p) = WSHED;
            end

        elseif (L(p) == MASK)
            L(p) = WSHED;

        end
    end

% Check for other minima at the same intensity level (i.e. more than one pixel has the same level: -> k2-k1~=1)

    for k = k1:k2
        p = ind(k);
        imd(p) = 0; % there are no neighbors, so set distance to 0

        if (L(p) == MASK) % then, p is inside a new minimum.
            currentLabel = currentLabel+1; % creating a new label
            queue(auxqueue) = p;
            auxqueue = auxqueue+1;
            L(p) = currentLabel;

            while ~isempty(queue)
                q = queue(1);
                queue(1) = [];
                auxqueue = auxqueue-1;

                s=NGneighbors(imi,q);

                for i=1:length(s);
                    r=s(i);

                    if (L(r) == MASK)
                        queue(auxqueue) = r;
                        auxqueue = auxqueue+1;
                        L(r) = currentLabel;
                    end

                end
            end
        end
    end
end
seg_im = L;


function q = NGneighbors(imi,p);
%
% function q = NGneighbors(imi,p);
%
% this function computes the neighbors pixels of p with 8-connectivity and
% return their linear (see help of sub2ind for further information) indexes
% into a single array q
%
% this function is called by: "water"
%
% author: Taian Vieira
%   date: 05/09/08


[x,y] = ind2sub(size(imi),p);
[y,x] = meshgrid([y-1 y y+1],[x-1 x x+1]);

% check if the pixel p is at any corner of the matrix
if numel([find(x==0 | x==size(imi,1)+1);find(y==0 | y==size(imi,2)+1)])==6
    x(x==0|x==size(imi,1)+1)=[]; x([1 2])=[];
    y(y==0|y==size(imi,2)+1)=[]; y([3 4])=[];
% otherwise the pixel may be along the first or last row
elseif ~isempty(find(x==0 | x==size(imi,1)+1))
    x(x==0|x==size(imi,1)+1)=[]; y(1,:)=[];
% or even along the first or last column
elseif ~isempty(find(y==0 | y==size(imi,2)+1))
    y(y==0|y==size(imi,2)+1)=[]; x(:,1)=[];
end

% converting the matrix index into a linear array of indexes
q = sub2ind(size(imi),x(:),y(:));
q(q==p)=[];


function equalized_im = equalization(im,RelativeEffortLevel,plotflag)

% the following two lines corresponds to a correction factor to define
% intensity extremities
MAXintensity = max(im(:))*RelativeEffortLevel;
im = round(im/MAXintensity*255); % creating integers of 8 bits according to the intensity full scale
F = cumsum(hist(im(:),[min(im(:)):max(im(:))]))/prod(size(im));
% F = cumsum(hist(im(:),linspace(0,MAXintensity,255)))/prod(size(im));

indexes = im - min(im(:)) + 1;
equalized_im = round([F(indexes) - min(F)] / [1 - min(F)] * 255 + 0.5);

% plotting comands
if nargin==2, return, end
f = cumsum(hist(im(:),0:255))/prod(size(im));
ax(1) = subplot(2,2,1),hist(im(:),0:255),hold,plot(0:255,f*max(hist(im(:),255)),'r')
subplot(2,2,3),imagesc(im)
f = cumsum(hist(equalized_im(:),0:255))/prod(size(im));
ax(2) = subplot(2,2,2),hist(equalized_im(:),0:255),hold,plot(0:255,f*max(hist(equalized_im(:),255)),'r')
subplot(2,2,4),imagesc(equalized_im)
set(ax,'xlim',[0 255])
% WATERSHED SET OF FUNCTIONS ------