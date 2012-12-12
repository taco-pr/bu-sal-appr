
close all;
clear all;


addpath('/home/robert/uni/matlabfns/NaN');
addpath('/home/robert/uni/taco/data/new');
addpath('/home/robert/uni/taco/svn.taco-project.eu/WP04/Software/SourceCode/Drafts/Libraries/matlab/NetCDF');


%figOrigRange = figure;
figOrigColor = figure;
figSaliency = figure;
subplotNM = [2 3];


%% configure saliency
p{1} = default_pami_param;
p{1}.mapWidth = 2;
p{1}.levels = [2 3 4];
p{1}.surroundSig = [2 8];
p{1}.channels = 'O';
p{1}.subtractMin = 1;
p{1}.useNormWeights = 0;
p{1}.nGaborAngles = 8;
p{1}.JumpEdgesThres = [0.8 0.9];
p{1}.RoofEdgesThres = [0.8 0.9];

% p{2} = p{1};
% p{2}.channels = 'R';
% 
% p{3} = p{1};
% p{3}.channels = 'S';
% 
% p{4} = p{1};
% p{4}.channels = 'T';
% 
% p{5} = p{1};
% p{5}.channels = 'U';
% 
% p{6} = p{1};
% p{6}.channels = 'V';

for nframe = 20:2:50

    %[img,NoRangeIndex] = getImageWithRangeFromCDF('kinect_with_RGB.cdf',nframe,1);
    [img,NoRangeIndex] = getImageWithRangeFromCDF('seq_003.cdf',nframe,1);

    figure(figOrigColor);imagesc(img(:,:,1:3));
    %figure(figOrigRange);imagesc(img(:,:,4));colormap(gray);

 
    for nparam = 1:numel(p)
        % calculate saliency
        [map{nparam} chanmaps{nparam} maps{nparam} chans{nparam}]= my_simpsal(img,p{nparam});
    end
    
    figure(figSaliency);clf;
    for nparam = 1:numel(p)
        % plot original image with overlayed saliency
        subplot(subplotNM(1),subplotNM(2),nparam);
        imagesc(img(:,:,1:3));hold on;
        
        if true
            h=imagesc(imresize(map{nparam},[size(img,1) size(img,2)],'nearest'));
            alpha_matrix = 0.9*ones(size(img,1), size(img,2));
            set(h,'AlphaData',alpha_matrix);
        else
            h=imagesc(0*ones(size(img,1), size(img,2)));       
            alpha_matrix = imresize(map{nparam},[size(img,1) size(img,2)],'nearest');
            alpha_matrix(alpha_matrix > 0.5) = 0;
            alpha_matrix(alpha_matrix ~= 0) = 0.8;
            set(h,'AlphaData',alpha_matrix);
        end
        
    end
end