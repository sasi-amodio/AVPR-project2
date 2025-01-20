
function image_with_nipples = nipple_amodio(image)

    % --------------------------------- gray scale image ---------------------------------- %
    
    bw_image = rgb2gray(image);
    
    %imshow(bw_image);
    
    % ------------------------------ human-body segmentation ------------------------------ %
    
    % human body binary mask 
    T = 50;
    [m,n] = size(bw_image);
    binary_mask = zeros([m,n] ,'logical');
    
    for i = 1 : m 
        for j = 1 : n
            if image(i,j) <= T
                binary_mask(i,j) = 0;
            else 
                binary_mask(i,j) = 1;
            end
        end
    end 
    
    % morphological closing
    SE = strel('disk',3);
    binary_mask =imclose(binary_mask, SE);
    
    % morphological dilation
    SE = strel('disk',10);
    human_body_binary_mask = imdilate(binary_mask,SE);
    
    %imshow(human_body_binary_mask);
    
    % -------------------------------- adaptive threshold -------------------------------- %
        
    C = 0.03;
    neighborhood = 15;
    mf_image  = medfilt2(bw_image,[neighborhood neighborhood]);
    image_tmp =  mf_image - bw_image;
    segmented_image = imbinarize(image_tmp,C);
    
    segmented_image = immultiply(segmented_image, human_body_binary_mask);
    
    %imshow(segmented_image);
    
    % -------------------------------------- mask_UL ------------------------------------- %

    [Him,Lim] = size(segmented_image);
    Lcnt = Lim / 2;
    Hup = 0.35 * Him;
    Hlw = 0.3 * Him;
    
    mask_UL = segmented_image(:,:);
     
    mask_UL(1:Hup,:) = 0;
    mask_UL((Him-Hlw):Him,:) = 0;
    
    %imshow(mask_UL);
    
    % --------------------------------- region detection --------------------------------- %
    
    
    %list of candidate computed only considerating valid regions 
    Np = 20;
    list_of_candidates = regionprops(mask_UL, 'PixelList', 'Circularity', 'Centroid');
    
    for i=1:size(list_of_candidates,1)
        candidate_size = size(list_of_candidates(i).PixelList,1);
        if (candidate_size < Np)
            list_of_candidates(i).PixelList=[]; 
            list_of_candidates(i).Centroid=[];
            list_of_candidates(i).Circularity=[];
        end
    end
    
    %regions_image = bw_image(:,:);
    %regions_image(:,:) = 0;
    %imshow(regions_image);
    %regions = cat(1,list_of_candidates.PixelList);
    %hold on
    %plot(regions(:,1),regions(:,2),'b*');
    %hold off
    
    
    % ---------------------------- nipple selection algorithm ---------------------------- %
     
    % distributes the candidates into L_left and L_right
    
    L_left = list_of_candidates;
    L_left_size = 0;
    
    L_right = list_of_candidates(:,:);
    L_right_size = 0;
    
    for i=1:size(list_of_candidates,1)
        candidate_centroid = list_of_candidates(i).Centroid;
        if(~isempty(candidate_centroid))
            if (candidate_centroid(1)> Lcnt) 
                L_left(i).PixelList=[]; 
                L_left(i).Centroid=[];
                L_left(i).Circularity=[];
                L_right_size = L_right_size+1;
            else 
                L_right(i).PixelList=[]; 
                L_right(i).Centroid=[];
                L_right(i).Circularity=[];   
                L_left_size = L_left_size+1;
            end
        end
    end
    
    % nipple selection algorithm

    if(L_left_size > 1) 
        % list of all the roundnesses of the candidates 
        candidates_roundnesses = [];
        % list of indices associated with candidates
        candidates_indeces = [];
    
        for i = 1:size(L_left,1)
            candidate_roundness = L_left(i).Circularity;
            if(~isempty(candidate_roundness))
                candidates_roundnesses = [candidates_roundnesses candidate_roundness];
                candidates_indeces = [candidates_indeces i];
            end
        end
    
        [~, m_l]= max(candidates_roundnesses);
        if(size(m_l) > 1)
            a_l = -1;
            for i=1:size(m_l) 
                % takes the index of each candidate associated with maximum roundness from the candidates_indeces
                index = candidates_indeces(m_l(i));
                candidate_area = L_left((index)).PixelList;
                if(candidate_area > a_l)
                    a_l = index;
                end
            end
            nipple_left = L_left(a_l);
        else
            index = candidates_indeces(m_l);
            nipple_left = L_left(index);
        end
    
    elseif (L_left_size == 1)
        for i=1:size(L_left,1)
            % the camndidate will be one
            candidate_centroid = L_left(i).Centroid;
            if(~isempty(candidate_centroid))
                nipple_left = L_left(i);
            end
        end
    else 
        nipple_left = [];
    end
    
    
    if(L_right_size > 1) 
        candidates_roundnesses = [];
        candidates_indeces = [];
        for i=1:size(L_right,1)
            candidate_roundness = L_right(i).Circularity;
            if(~isempty(candidate_roundness))
                candidates_roundnesses = [candidates_roundnesses candidate_roundness];
                candidates_indeces = [candidates_indeces i];
            end
        end
    
        [~, m_r] = max(candidates_roundnesses);
        if(size(m_r) > 1)
            a_r = 0;
            for i=1:size(m_r) 
                % takes the index of each candidate associated with maximum roundness from the candidates_indeces
                index = candidates_indeces(m_r(i));
                candidate_area = L_right((index)).PixelList;
                if(candidate_area > a_r)
                    a_r = index;
                end
            end
            nipple_right = L_right(a_r);
        else
            index = candidates_indeces(m_r);
            nipple_right = L_right(index);
        end
    elseif (L_right_size == 1)
        for i=1:size(L_right,1)
            % the candidate will be one
            candidate_centroid = L_right(i).Centroid;
            if(~isempty(candidate_centroid))
                nipple_right = L_right(i);
            end
        end
    else 
        nipple_right = [];
    end
    

    image_with_nipples = bw_image(:,:);
    imshow(image_with_nipples);

    hold on
    if(~isempty(nipple_left))
        plot(nipple_left.Centroid(:,1),nipple_left.Centroid(:,2),'b*', 'LineWidth',1.2, 'Markersize', 7);
    end
    if(~isempty(nipple_right))
        plot(nipple_right.Centroid(:,1),nipple_right.Centroid(:,2),'b*', 'LineWidth',1.2, 'Markersize', 7);
    end
    print('-dpng','-r0', 'results.png');
    hold off    
    
    image_with_nipples = imread('results.png');
end