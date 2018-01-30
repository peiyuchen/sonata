function [ S_matrix ] = similar_matrix( Pattern1, Pattern2, distance, S_plot )
% ��ǭ�bar�Z��,�o��similarity matrix
% input  : �����pattern, �Ϊ�distance
% output : �����x�}

    if nargin < 3, distance = 'correlation'; end
    if nargin < 4, S_plot = 0;               end
    
    %% distance
    D_matrix    = pdist2(Pattern1', Pattern2', distance); 
%     twoVecDist  = pdist2(Pattern', distance);        % ��Ǥ������Z��
%     D_matrix    = squareform(twoVecDist);           % �ܦ�distance Matrix
    
    %% similarity
    S_matrix = 1 - D_matrix;
    S_matrix(S_matrix < 0) = 0;

    %% plot similarity matrix
    if S_plot
        figure, imagesc(S_matrix); colorbar;
        title(['similar , distance = ' distance]);
    end    
    

end

