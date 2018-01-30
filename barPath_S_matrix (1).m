function [ S_matrix ] = barPath_S_matrix( Pattern1, Pattern2, distance, S_plot )
% 算倆倆bar距離,得到similarity matrix,    多音bar 會有多個路徑
% input  : 比較的pattern, 用的distance
% output : 相關矩陣

    if nargin < 3, distance = 'correlation'; end
    if nargin < 4, S_plot = 0;               end
    
    h = waitbar(0,'Please wait...');
          
    for i=1:length(Pattern1.barPath)
        barPath1 = Pattern1.barPath{i,:};
        
        for j =1:length(Pattern2.barPath)
            barPath2 = Pattern2.barPath{j,:};
%             [1,i,j]
            
            tmp_S_matrix = similar_matrix(barPath1, barPath2, distance);
%             [2,i,j]
            
            max1 = max(tmp_S_matrix,[],1);
            max2 = max(tmp_S_matrix,[],2);

            
            S{i,j}        = tmp_S_matrix;
            S_matrix(i,j) = mean([max1 max2']);
            

%             if i==9 && j == 10
%                 figure, 
%                 subplot(2,2,1); plot(barPath1); axis([0 inf 50 70]); title(['Path : bar' num2str(i)]); legend('show'); xlabel('beat'); ylabel('pitch');
%                 subplot(2,2,2); plot(barPath2); axis([0 inf 50 70]); title(['Path : bar' num2str(j)]); legend('show'); xlabel('beat'); ylabel('pitch');
%                 subplot(2,2,3:4);
%                 imagesc(tmp_S_matrix);  title(['Similar\_matrix :[' num2str(i) ',' num2str(j) ']']);  colormap(flipud(gray)); colorbar;
%             end  
%             if i==9 && j == 12
%                 figure, 
%                 subplot(2,2,1); plot(barPath1); axis([0 inf 50 70]); title(['Path : bar' num2str(i)]); legend('show'); xlabel('beat'); ylabel('pitch');
%                 subplot(2,2,2); plot(barPath2); axis([0 inf 50 70]); title(['Path : bar' num2str(j)]); legend('show'); xlabel('beat'); ylabel('pitch');
%                 subplot(2,2,3:4);
%                 imagesc(tmp_S_matrix);  title(['Similar\_matrix']);  colormap(flipud(gray)); colorbar;
%             end 
            
        end
        waitbar(i/length(Pattern1.barPath));
    end
    close(h);
    
    %% plot similarity matrix
    if S_plot
        figure, imagesc(S_matrix); colorbar;
        title(['similar , distance = ' distance]);
    end  

end

    
%     
%     if plotPattern
%         figure, imagesc(S_matrix);  title(['S\_matrix\_multiPitch']);  colormap(flipud(gray)); 
%     end  

