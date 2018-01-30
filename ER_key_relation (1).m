function [ PATHFAMILY,SHIFTMATRIX ] = ER_key_relation( Path1, Path2, pathFamily, filename, multiProcess, p, S_no )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
    [GT_subject] = structure_annotation(['../data_annotation/pei_anno/' filename '.txt'], 2);
    
    f1 = figure;
    PATHFAMILY.SM               = similar_matrix(Path1, Path2, 'corr');
    SHIFTMATRIX.shift_matrix    = zeros(12,size(Path1, 2));

    for i=1:size(Path1, 2)
        % 大於0.85相似度的才去看
        if similar_matrix(Path1(:,i), Path2(:,i), 'corr') > 0.85
            PATHFAMILY.shift(:,i) = mode(mod(round(Path1(:,i) - Path2(:,i)),12));
%             if i>49&&i<=49*2
% %                 if i==1;figure;end
%                 subplot(7,7,i-49);
%                 plot(Path1(:,i)); hold all; plot(Path2(:,i)); title([num2str(pathFamily(1,i)) '&' num2str(pathFamily(2,i))]);
%             end

            % 畫 shift 圖
            plot(pathFamily(1,i), PATHFAMILY.shift(:,i),'.','MarkerSize',20,'Color','b'); hold all; GT_line(GT_subject, 2*pathFamily(1,1));
            xlabel('Time (bar)');  ylabel('Shift index');  axis([pathFamily(1,1), pathFamily(1,end), 0, 12]);  title(['no.' filename(3:4)]);

            % 畫成 12*N shift matrix
            SHIFTMATRIX.shift_matrix(PATHFAMILY.shift(:,i)+1,i) = 1;
        end
        
    end
    
    %% plot
    % Erode & Dilate
%     span    = 5;
%     SE      = strel('line',span,0);
%     Dilate  = imdilate(SHIFTMATRIX.shift_matrix,SE);
%     Erode   = imerode (SHIFTMATRIX.shift_matrix,SE);
%     f3      = figure;
%     subplot(311);   plotShiftMatrix(SHIFTMATRIX.shift_matrix, 'Normalize', pathFamily, GT_subject );    
%     subplot(312);   plotShiftMatrix(                  Dilate,    'Dilate', pathFamily, GT_subject );
%     subplot(313);   plotShiftMatrix(                   Erode,     'Erode', pathFamily, GT_subject );
%     f4      = figure;
%     subplot(311);   plotShift(SHIFTMATRIX.shift_matrix, 'Normalize', pathFamily, GT_subject );    
%     subplot(312);   plotShift(                  Dilate,    'Dilate', pathFamily, GT_subject );
%     subplot(313);   plotShift(                   Erode,     'Erode', pathFamily, GT_subject );


   % closing opening
    tmp       = SHIFTMATRIX.shift_matrix;    
    span      = 3;
    SE        = strel('line',span,0);
    Opening   = imdilate(imerode(tmp,SE),SE);
    Closing   = imerode(imdilate(tmp,SE),SE);
    OC        = imerode(imdilate(Opening,SE),SE);
    CO        = imdilate(imerode(Closing,SE),SE);
    Morphology= figure;
    subplot(5,1,1); plotShiftMatrix(    tmp, 'Normalize'                                              , pathFamily, GT_subject );
    subplot(5,1,2); plotShiftMatrix(Opening, ['Erode -> Dilate, ' num2str(length(SE.Neighborhood))]   , pathFamily, GT_subject );          
    subplot(5,1,3); plotShiftMatrix(Closing, ['Dilate -> Erode, ' num2str(length(SE.Neighborhood))]   , pathFamily, GT_subject );
    subplot(5,1,4); plotShiftMatrix(     OC, ['Opening -> Closing, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
    subplot(5,1,5); plotShiftMatrix(     CO, ['Closing -> Opening, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
    
    Morphology1= figure;
    subplot(5,1,1); plotShift(    tmp, 'Normalize'                                              , pathFamily, GT_subject );
    subplot(5,1,2); plotShift(Opening, ['Erode -> Dilate, ' num2str(length(SE.Neighborhood))]   , pathFamily, GT_subject );          
    subplot(5,1,3); plotShift(Closing, ['Dilate -> Erode, ' num2str(length(SE.Neighborhood))]   , pathFamily, GT_subject );
    subplot(5,1,4); plotShift(     OC, ['Opening -> Closing, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
    subplot(5,1,5); plotShift(     CO, ['Closing -> Opening, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
    
    f5 = figure;
    subplot(2,1,1); plotShift(    tmp, 'Normalize'                                              , pathFamily, GT_subject );
    subplot(2,1,2); plotShift(     OC, ['Opening -> Closing, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
%     saveas(        f5, ['ER_relation/' multiProcess '/no.' filename(3:4) '_OC_shift_Sno_' num2str(S_no) '_pf' num2str(p) '.png']);
%     
%     saveas(        f1 , ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift_Sno' num2str(S_no) '_pf' num2str(p) '.png']);
%     saveas(Morphology , ['ER_relation/' multiProcess '/no.' filename(3:4) '_Morphology_Matrix_Sno' num2str(S_no) '_pf' num2str(p) '.png']);
%     saveas(Morphology1, ['ER_relation/' multiProcess '/no.' filename(3:4) '_Morphology_Sno' num2str(S_no) '_pf' num2str(p) '.png']);
    saveas(f5, ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift' num2str(S_no) '_pf' num2str(p) '.png']);
    
end
 

function plotShift(data, name, pathFamily, GT_subject)
        for i=1:size(data, 2)
            idx         = find(data(:,i)==1)-1;
            if ~isempty(idx)
                x = repmat(pathFamily(1,i),1,length(idx));
                plot(x, idx,'.','MarkerSize',20,'Color','b'); hold all;
            end
        end
        xlabel('Time (bar)');  ylabel('Shift index');  axis([pathFamily(1,1), pathFamily(1,end), 0, 12]);  title(name);      
        GT_line(GT_subject, pathFamily);
end

function GT_line(GT_subject, pathFamily)
    for GT = 1:length(GT_subject)
        line([GT_subject(GT).M1_start GT_subject(GT).M1_start],[0 12],'Color','g','LineWidth',1); hold all;
        line([GT_subject(GT).M1_end   GT_subject(GT).M1_end]  ,[0 12],'Color','g','LineWidth',2); hold all;
        line([GT_subject(GT).M2_start GT_subject(GT).M2_start],[0 12],'Color','r','LineWidth',2); hold all;
        line([GT_subject(GT).M2_end   GT_subject(GT).M2_end]  ,[0 12],'Color','r','LineWidth',1); hold all;
    end
end

function plotShiftMatrix(data, name, pathFamily, GT_subject)
        imagesc(abs(data-1)); title(name); hold all; colormap gray;
        set(gca,'YDir','normal'); set(gca,'YTick',[1:2:12]); set(gca,'YTickLabel',0:2:11);
        set(gca,'XTick', [1:10:length(pathFamily)]);  set(gca,'XTickLabel', [pathFamily(1,1):10:pathFamily(1,end)]); 
        GT_line1(GT_subject, pathFamily);
end

function GT_line1(GT_subject, pathFamily)
    for GT = 1:length(GT_subject)
        if GT_subject(GT).M1_start-pathFamily(1,1)+2>0
            line([GT_subject(GT).M1_start-pathFamily(1,1)+1 GT_subject(GT).M1_start-pathFamily(1,1)+1],[0 12],'Color','g','LineWidth',1); hold all;
        end
        if GT_subject(GT).M1_end-pathFamily(1,1)+1>0
            line([GT_subject(GT).M1_end-pathFamily(1,1)+1   GT_subject(GT).M1_end-pathFamily(1,1)+1]  ,[0 12],'Color','g','LineWidth',2); hold all;
            line([GT_subject(GT).M2_start-pathFamily(1,1)+1 GT_subject(GT).M2_start-pathFamily(1,1)+1],[0 12],'Color','r','LineWidth',2); hold all;
            line([GT_subject(GT).M2_end-pathFamily(1,1)+1   GT_subject(GT).M2_end-pathFamily(1,1)+1]  ,[0 12],'Color','r','LineWidth',1); hold all;
        end
    end
end