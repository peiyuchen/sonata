function [ S_matrix ] = similar_matrix( Pattern1, Pattern2, distance, S_plot )
% 算倆倆bar距離,得到similarity matrix
% input  : 比較的pattern, 用的distance
% output : 相關矩陣

    if nargin < 3, distance = 'correlation'; end
    if nargin < 4, S_plot = 0;               end
    
    %% distance
    D_matrix    = pdist2(Pattern1', Pattern2', distance); 
%     twoVecDist  = pdist2(Pattern', distance);        % 兩倆之間的距離
%     D_matrix    = squareform(twoVecDist);           % 變成distance Matrix
    
    %% similarity
    S_matrix = 1 - D_matrix;
    S_matrix(S_matrix < 0) = 0;

    %% plot similarity matrix
    if S_plot
        figure, imagesc(S_matrix); colorbar;
        title(['similar , distance = ' distance]);
    end    
    

end

