function [ PATHFAMILY,SHIFTMATRIX ] = ER_key_relation_v2( Path1, Path2, pathFamily, filename, multiProcess )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
%     [GT_subject] = parseAnnotationMovementFile(['../data_annotation/' filename '_movement.txt']); % 讀 主題的標記  

    [GT_subject] = structure_annotation(['../data_annotation/pei_anno/' filename '.txt'], 2);
    
%     PATHFAMILY.shift_matrix = zeros(12,size(Path1, 1));
    
    f1 = figure;
    SHIFTMATRIX.shift_matrix = zeros(12,size(Path1, 1));

    for i=1:size(Path1, 1)
        
        % 對應的bar 所有路徑比較相似性
        PATHFAMILY(i).SM = similar_matrix(Path1{i,:}, Path2{i,:}, 'corr');
         
        %{ 
%       路徑較少的bar 找與這個bar的所有路徑最像的bar(取max)，且相似度要大於0.85
%         if size(PATHFAMILY.SM{i,:}, 1) <= size(PATHFAMILY.SM{i,:}, 2)
%             [PATHFAMILY.SM_MAX1{i}, id] = max(PATHFAMILY.SM{i,:}, [], 2);
%             delete = find(PATHFAMILY.SM_MAX1{i}<0.85)';
% 
%             PATHFAMILY.MAX1_p{i,1}   = Path1{i,:};
%             PATHFAMILY.MAX1_p{i,2}   = Path2{i}(:,id);  
%             PATHFAMILY.MAX1_p{i,1}(:,delete) = [];
%             PATHFAMILY.MAX1_p{i,2}(:,delete) = [];
%         else
%             [PATHFAMILY.SM_MAX1{i}, id] = max(PATHFAMILY.SM{i,:}, [], 1);
%             delete = find(PATHFAMILY.SM_MAX1{i}<0.85)';
% 
%             PATHFAMILY.MAX1_p{i,1}   = Path1{i}(:,id);
%             PATHFAMILY.MAX1_p{i,2}   = Path2{i,:};
%             PATHFAMILY.MAX1_p{i,1}(:,delete) = [];
%             PATHFAMILY.MAX1_p{i,2}(:,delete) = [];
%         end
%}
        % 取max，相似度要大於0.85
            [PATHFAMILY(i).SM_MAX1, id] = max(PATHFAMILY(i).SM, [], 2);
            delete = find(PATHFAMILY(i).SM_MAX1<0.85)';
            PATHFAMILY(i).delete = delete;
            
            PATHFAMILY(i).MAX1_p{1}   = Path1{i,:};
            PATHFAMILY(i).MAX1_p{2}   = Path2{i}(:,id);  
            PATHFAMILY(i).MAX1_p{1}(:,delete) = [];
            PATHFAMILY(i).MAX1_p{2}(:,delete) = [];
            
            [PATHFAMILY(i).SM_MAX2, id2] = max(PATHFAMILY(i).SM, [], 1);
            delete2 = find(PATHFAMILY(i).SM_MAX2<0.85)';
            PATHFAMILY(i).delete2 = delete2;
            PATHFAMILY(i).MAX2_p{1}   = Path1{i}(:,id2);
            PATHFAMILY(i).MAX2_p{2}   = Path2{i,:};
            PATHFAMILY(i).MAX2_p{1}(:,delete2) = [];
            PATHFAMILY(i).MAX2_p{2}(:,delete2) = [];

        % 留下了的路徑算shift
        PATHFAMILY(i).shift1 = [];
        PATHFAMILY(i).shift2 = [];

        PATHFAMILY(i).shift1      = mode(mod(round(PATHFAMILY(i).MAX1_p{1} - PATHFAMILY(i).MAX1_p{2}+12), 12));
        PATHFAMILY(i).shift2      = mode(mod(round(PATHFAMILY(i).MAX2_p{1} - PATHFAMILY(i).MAX2_p{2}+12), 12));
        if isnan(PATHFAMILY(i).shift1); PATHFAMILY(i).shift1 = []; end
        if isnan(PATHFAMILY(i).shift2); PATHFAMILY(i).shift2 = []; end
        PATHFAMILY(i).shift       = [PATHFAMILY(i).shift1 PATHFAMILY(i).shift2];

        % shift index的個數
        shift_element       = unique(PATHFAMILY(i).shift + 1);
        count_shift_element = hist(PATHFAMILY(i).shift+1, shift_element);
        if length(shift_element)==1; count_shift_element = sum(count_shift_element); end

        % 畫 shift 圖
        if ~isempty(PATHFAMILY(i).shift)
            plot(pathFamily(1,i), PATHFAMILY(i).shift,'.','MarkerSize',20,'Color','b'); hold all;
            xlabel('Time (bar)');   ylabel('Shift index');  axis([pathFamily(1,1), pathFamily(1,end), 0, 12]);  title(['no.' filename(3:4)]);
            GT_line(GT_subject, pathFamily);
        end
        % 畫成 12*N shift matrix
        SHIFTMATRIX.shift_matrix((i-1)*12 + shift_element) = count_shift_element + SHIFTMATRIX.shift_matrix((i-1)*12 + shift_element);

%         PATHFAMILY(i).shift_matrix( (i-1)*12 + shift_element) = count_shift_element + PATHFAMILY(i).shift_matrix( (i-1)*12 + shift_element);
        if max(SHIFTMATRIX.shift_matrix(:,i)) ~= 0
            SHIFTMATRIX.shift_matrix_nor(:,i)               = SHIFTMATRIX.shift_matrix(:,i)./max(SHIFTMATRIX.shift_matrix(:,i));
        end    
    end
    
   
    % shift matrix 做 smooth
    span = round(length(Path1)/10);
    for i =1:12
        SHIFTMATRIX.shift_smooth(i,:) = smooth(SHIFTMATRIX.shift_matrix_nor(i,:), span); 
    end

    % threshold
    th = 0.5;
    SHIFTMATRIX.shift_th = SHIFTMATRIX.shift_smooth;
    SHIFTMATRIX.shift_th( SHIFTMATRIX.shift_th < th ) = 0;
    
    f2 = figure; 
    subplot(221);   plotShiftMatrix(SHIFTMATRIX.shift_matrix, 'shift\_matrix', pathFamily, GT_subject );  
    subplot(222);   plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'normalize', pathFamily, GT_subject );  
    subplot(223);   plotShiftMatrix(SHIFTMATRIX.shift_smooth, ['smooth, span = ' num2str(span)], pathFamily, GT_subject );  
    subplot(224);   plotShiftMatrix(SHIFTMATRIX.shift_th, ['threshold = ' num2str(th)], pathFamily, GT_subject );  

        % Erode & Dilate
        span    = 5;
        SE      = strel('line',span,0);
        Dilate  = imdilate(SHIFTMATRIX.shift_matrix_nor,SE);
        Erode   = imerode(SHIFTMATRIX.shift_matrix_nor,SE);
        f3      = figure;
        subplot(311);   plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );           
        subplot(312);   plotShiftMatrix(Dilate, 'Dilate', pathFamily, GT_subject );    
        subplot(313);   plotShiftMatrix(Erode, 'Erode', pathFamily, GT_subject );    

        DE      = imerode(Dilate,SE);
        ED      = imdilate(Erode,SE);
        f4      = figure;
        subplot(311);   plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );    
        subplot(312);   plotShiftMatrix(DE, 'Dilate->Erode', pathFamily, GT_subject );     
        subplot(313);   plotShiftMatrix(ED, 'Erode->Dilate', pathFamily, GT_subject );             
        
        tmp = SHIFTMATRIX.shift_matrix_nor;
        figure;
        subplot(6,1,1); plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );             
        for i = 2:6
            span = 2*i-1;
            SE = strel('line',span,0);
            TTT = imerode(imdilate(tmp,SE),SE);
            subplot(6,1,i); plotShiftMatrix(TTT, ['Dilate -> Erode, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );             
        end
        
        tmp = SHIFTMATRIX.shift_matrix_nor;
        figure;
        subplot(6,1,1); plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject ); 
        for i = 2:6
            span = 2*i-1;
            SE = strel('line',span,0);
            TTT = imdilate(imerode(tmp,SE),SE);
            subplot(6,1,i); plotShiftMatrix(TTT, ['Erode -> Dilate, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject ); 
            
        end
        
        
        tmp = SHIFTMATRIX.shift_matrix_nor;
        figure;
        subplot(6,1,1); plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );        
        for i = 2:6
            span = 3;
            SE = strel('line',span,0);
            TTT = imdilate(imerode(tmp,SE),SE);
            tmp = TTT;
            subplot(6,1,i); plotShiftMatrix(TTT, ['Erode -> Dilate, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );  
        end
        
        tmp = SHIFTMATRIX.shift_matrix_nor;
        figure;
        subplot(6,1,1);   plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );
        for i = 2:6
            span = 2*i-1;
            SE = strel('line',span,0);
            TTT = imerode(imdilate(imdilate(imerode(tmp,SE),SE),SE),SE);
%             tmp = TTT;
            subplot(6,1,i); plotShiftMatrix(TTT, ['E -> D -> D -> E, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );           
        end
        
        
        tmp = SHIFTMATRIX.shift_matrix_nor;
        figure;
        subplot(6,1,1);   plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );  
        for i = 2:6
            span = 2*i-1;
            SE = strel('line',span,0);
            TTT = imdilate(imerode(imerode(imdilate(tmp,SE),SE),SE),SE);
            subplot(6,1,i); plotShiftMatrix(TTT, ['D -> E -> E -> D, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );  
            
        end
        
        
        
       tmp =SHIFTMATRIX.shift_matrix_nor; T_NAME = []; Y_NAME = []; tmp1 = tmp;
        DE_f = figure;
        subplot(7,2,1:2);  plotShiftMatrix(SHIFTMATRIX.shift_matrix_nor, 'Normalize', pathFamily, GT_subject );  
        for i = 2:7
            if (i==2 || i==3);    span = 3;   end
            if (i==4 || i==5);    span = 5;   end
            if mod(i,2) == 0
                span = i+1;
            else
                span = i;
            end
            SE = strel('line',span,0);
            if mod(i,2) == 0
                TTT = imerode(tmp,SE);
                YYY = imdilate(tmp1,SE);
                t_name = 'E';
                y_name = 'D';
            else
                TTT = imdilate(tmp,SE);
                YYY = imerode(tmp1,SE);
                t_name = 'D';
                y_name = 'E';
                
            end

                tmp = TTT; if i == 2; T_NAME = t_name; else T_NAME = [T_NAME '->' t_name]; end
                subplot(7,2,2*i-1);   plotShiftMatrix(TTT, [num2str(length(SE.Neighborhood)) T_NAME], pathFamily, GT_subject );
                tmp1 = YYY; if i == 2; Y_NAME = y_name; else Y_NAME = [Y_NAME '->' y_name]; end
                subplot(7,2,2*i);     plotShiftMatrix(YYY, [num2str(span) Y_NAME], pathFamily, GT_subject );    
        end
        
        
        
        tmp       = SHIFTMATRIX.shift_matrix_nor;    
        span      = 3;
        SE        = strel('line',span,0);
        Opening   = imdilate(imerode(tmp,SE),SE);
        Closing   = imerode(imdilate(tmp,SE),SE);
        OC        = imerode(imdilate(Opening,SE),SE);
        CO        = imdilate(imerode(Closing,SE),SE);
        Morphology= figure;
        subplot(5,1,1); plotShiftMatrix(tmp, 'Normalize', pathFamily, GT_subject );
        subplot(5,1,2); plotShiftMatrix(Opening, ['Erode -> Dilate, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );          
        subplot(5,1,3); plotShiftMatrix(Closing, ['Dilate -> Erode, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
        subplot(5,1,4); plotShiftMatrix(OC, ['Opening -> Closing, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
        subplot(5,1,5); plotShiftMatrix(CO, ['Closing -> Opening, ' num2str(length(SE.Neighborhood))], pathFamily, GT_subject );
        
%         saveas(f1, ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift.png']);
%         saveas(f2, ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift_matrix_smooth.png']);
%         saveas(f3, ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift_E&D.png']); 
%         saveas(f4, ['ER_relation/' multiProcess '/no.' filename(3:4) '_shift_ED&DE.png']); 
%         saveas(DE_f, ['ER_relation/' multiProcess '/no.' filename(3:4) 'iter_DE.png']); 
%         saveas(Morphology, ['ER_relation/' multiProcess '/no.' filename(3:4) '_Morphology_span3.png']); 


end

function GT_line(GT_subject, pathFamily)
    for GT = 1:length(GT_subject)
        if GT_subject(GT).M1_start-pathFamily(1,1)+1>0
            line([GT_subject(GT).M1_start-pathFamily(1,1)+1 GT_subject(GT).M1_start-pathFamily(1,1)+1],[0 12],'Color','g','LineWidth',1); hold all;
        end
        if GT_subject(GT).M1_end-pathFamily(1,1)+1>0
            line([GT_subject(GT).M1_end-pathFamily(1,1)+1 GT_subject(GT).M1_end-pathFamily(1,1)+1],[0 12],'Color','g','LineWidth',2); hold all;
            line([GT_subject(GT).M2_start-pathFamily(1,1)+1 GT_subject(GT).M2_start-pathFamily(1,1)+1],[0 12],'Color','r','LineWidth',2); hold all;
            line([GT_subject(GT).M2_end-pathFamily(1,1)+1 GT_subject(GT).M2_end-pathFamily(1,1)+1],[0 12],'Color','r','LineWidth',1); hold all;
        end
    end
end

function plotShiftMatrix(data, name, pathFamily, GT_subject)
        imagesc(abs(data-1)); title(name); hold all; colormap gray;
        set(gca,'YDir','normal'); set(gca,'YTick',[1:2:12]); set(gca,'YTickLabel',0:2:11);
        set(gca,'XTick', [1:10:length(pathFamily)]);  set(gca,'XTickLabel', [pathFamily(1,1):10:pathFamily(1,end)]); 
        GT_line(GT_subject, pathFamily);
end
% % p=1;
% % for i=1:length( pathFamily{p,1} ) 
% %     SS(i) = S_final( pathFamily{p,1}(2,i), pathFamily{p,1}(1,i) );
% % end
% % TH = PATHFAMILY(i)(1).both_shift_th;
% % TH(TH>0)=1;
% % 
% % SS.*TH;
