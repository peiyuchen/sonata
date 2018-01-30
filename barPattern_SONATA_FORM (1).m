%% ANALYZING MUSIC RECORDINGS IN SONATA FORM

% case 1 : AUTOMATED METHODS FOR ANALYZING MUSIC RECORDINGS IN SONATA FORM (PAPER 2013) 
%          MIDI FILE -> CENS feature -> Self-similarity matrix -> Thumbnail
% case 2 : combine the "bar pattern"
%          MIDI FILE -> Bar pattern  -> Self-similarity matrix -> Thumbnail


clear all; close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  0. Loads a midi file & Choose method : case 1 or case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = [1,2,3,4,5,8,16,17,20]; % 8


for song_no = 1:length(index)
clearvars -except index song_no; close all; clc;
filename = ['b_' num2str(index(song_no)) '_1']
% filename = ['b_' num2str(song_no) '_1'];

SAVE    =  0;
S_no    =  1;
multiProcess   = 'modifySkyline';  % new   , 'allPath','lowerPitch';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1. add paths of toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('toolbox/MATLAB_SM-Toolbox_1.0/SMtoolbox_functions');
% addpath('toolbox/MATLAB_SM-Toolbox_1.0/MATLAB-Chroma-Toolbox_2.0');
% addpath('toolbox/MIDI tool/miditoolbox');
% addpath('toolbox/midi_lib/midi_lib');
% addpath('toolbox/matlab-midi-master/src');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CASE 2 : Computes Bar pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(['S_matrix/' multiProcess '/' filename '.mat'])
    save_S_matrix( filename, multiProcess );
end
    load(['S_matrix/' multiProcess '/'  filename '.mat']);

    
    
S_name  = {'max3', 'max2', 'right', 'left', 'com'};

for S_no=3%1:5

    % 找 Smax index 
    if S_no == 1
        tmp_r = S.max3 - S.right; tmp_r(tmp_r==0)=-1; tmp_r(tmp_r~=-1)=0;
        tmp_l = S.max3 - S.left;  tmp_l(tmp_l==0)=2;  tmp_l(tmp_l~=2)=0;
        tmp_c = S.max3 - S.com;   tmp_c(tmp_c==0)=3;  tmp_c(tmp_c~=3)=0;
        S.index = zeros(size(S.max3)); S.index(tmp_c==3)=3; S.index(tmp_l==2)=2; S.index(tmp_r==-1)=1;
    elseif S_no == 2
        tmp_r = S.max3 - S.right; tmp_r(tmp_r==0)=-1; tmp_r(tmp_r~=-1)=0;
        tmp_l = S.max3 - S.left;  tmp_l(tmp_l==0)=2;  tmp_l(tmp_l~=2)=0;
        S.index = zeros(size(S.max2)); S.index(tmp_l==2)=2; S.index(tmp_r==-1)=1;
    elseif S_no == 3
        S.index = ones(size(S.right)); 
    elseif S_no == 4
        S.index = 2 * ones(size(S.left)); 
    elseif S_no == 5
        S.index = 3 * ones(size(S.com));
    end

    S_th = eval( [ 'S.' S_name{S_no}; ] );   
    % figure, imagesc(S_th);    title( S_name{S_no} );   % colormap(flipud(gray)); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   3. Computes and visualizes an thresholded similarity matrix. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    th                  = 0.85; % 取門檻值
    S_th( S_th < th )   = -2  ;
    S_th( S_th >= th )  = (S_th(S_th>=th)-min(min(S_th(S_th>=th)))) / ( max(max(S_th(S_th>=th))) - min(min(S_th(S_th>=th))) );
    % figure, imagesc(S_th);  title(['S\_th = ' num2str(th)]);  colormap(flipud(gray)); 

    S_final = S_th;



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   4. Computes and saves a fitness scape plot.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute fitness scape plot and save
    parameter.dirFitness    = ['data_fitness/' multiProcess '/'];
    parameter.saveFitness   = 1;
    parameter.title         = filename;
    parameter.fitFileName   = [filename '_S' S_name{S_no} '_fit'];
    parameter.fitnessSaveFileName = [parameter.dirFitness, filename,'_S' S_name{S_no} '_fit','.mat'];
    %fitnessSaveFileName = [parameter.dirFitness,parameter.fitFileName,'.mat'];

    %-----------!!IMPORTANT!!--------------------------------------------------%
    % For fast computing of fitness scape plot, please enable parallel computing.
    % To enable that, use command 'matlabpool open'.
    % To disable that, use command 'matlabpool close'
    %--------------------------------------------------------------------------%
    
    if ~exist([parameter.fitnessSaveFileName])        
        disp('NO fitness .mat');
        [fitness_info,parameter] = SSM_to_scapePlotFitness(S_final, parameter);
        fitness_matrix           = fitness_info.fitness;
    else % % instead of computing fitness, you can load a previously computed scape plot:
        disp('exist fitness .mat');
        fitnessFile = load(parameter.fitnessSaveFileName);
        fitness_matrix = fitnessFile.fitness_info.fitness;
    end

    paramVisScp = [];
    % paramVisScp.timeLineUnit = 'sample';
    % paramVisScp.timeLineUnit = 'second'; paramVisScp.featureRate = ... 
    [h_fig_scapeplot,x_axis,y_axis] = visualizeScapePlot(fitness_matrix,paramVisScp);
    title('Fitness scape plot','Interpreter','none');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  5. Computes the thumbnail 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if S_no ~= 5
    % compute thumbnail with length constraint
     parameter.len_min_seg_frame = ceil(size(S_final,1)/6);   % peiyu
    % parameter.len_min_seg_frame= 20;
    [thumb_frame] = scapePlotFitness_to_thumbnail(fitness_matrix,parameter);

    % show corresponding thumbnail point in fitness scape plot
    center_thumb_frame = floor((thumb_frame(1) + thumb_frame(2))/2);
    length_thumb_frame = thumb_frame(2) - thumb_frame(1) + 1;

    point_x_pos = x_axis(center_thumb_frame);
    point_y_pos = y_axis(length_thumb_frame);

    hold on;
    plot(point_x_pos,point_y_pos,'o','LineWidth',2,'color',[1, 0.5, 0]);
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   6. Computes optimal path family and induced segment family
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % find repetitions of thumbnail
    [induced_frame,pathFamily] = thumbnailSSM_to_pathFamily(thumb_frame,S_final,parameter);
    for p = 1:length(pathFamily)
        S_now = eval( [ 'S.' S_name{S_no} ]);
        D{p} = S_now(pathFamily{p,1}(1,1):pathFamily{p,1}(1,end), pathFamily{p,1}(2,1):pathFamily{p,1}(2,end));
        accu_M{p} = accumulate_matrix( D{p} );
        [d, ~] = dtw(accu_M{p});
        d(:,1) = d(:,1) + pathFamily{p,1}(1,1) - 1;
        d(:,2) = d(:,2) + pathFamily{p,1}(2,1) - 1;
        DTWpathFamily{p,1} = d';
    end
    
    

    paramVisPathSSM                          = [];
    paramVisPathSSM.visualizeInducedSegments =  1;
    paramVisPathSSM.visualizeWarpingpath     =  1;
    paramVisPathSSM.dirFigure                = ['data_image/'  multiProcess '/'];
    paramVisPathSSM.print                    = 1;
    paramVisPathSSM.figureName               = [filename '_pathFamily_BAR'];
    visualizePathFamilySSM(S_final,DTWpathFamily,paramVisPathSSM);
    title('S, path family, and induced segment family');


    % convert from frames to seconds
    %{
    % parameter.featureRate = 10/paramCENS.downsampSmooth;
    % parameter.duration = size(S_final,1)/parameter.featureRate;
    % induced_second = convertSegment_frames_to_seconds(induced_frame,parameter.featureRate);
    % thumb_second = convertSegment_frames_to_seconds(thumb_frame,parameter.featureRate);
    %}

    % attach audio file to SSMPathFamily
    %{
    % if isfield(parameter,'timeLineUnit') && (strcmp(parameter.timeLineUnit,'second'))
    %     parameterMPP.featureTimeResType = 'seconds';
    % else
    %     parameterMPP.featureTimeResType = 'features';
    % end
    % parameterMPP.featureRate = parameter.featureRate;
    % parameterMPP.fs = sideinfo.wav.fs;
    % h_fig = gcf;
    % makePlotPlayable(f_audio, h_fig, parameterMPP);
    %}

    % assign label to each repetition and wrap up in segment struct
    computedSegments = wrapUpSegmentInStruct(induced_frame,thumb_frame);   % peiyu
    % computedSegments = wrapUpSegmentInStruct(induced_second,thumb_second);
    paramVisSegFam = [];
    paramVisSegFam.duration = size(S_final,1);
    % paramVisSegFam.duration = parameter.duration;  % 畫圖x最大範圍

    % 畫計算過後的結果
    %{
    % paramVisSegFam.showLabelText = 1;
    % paramVisSegFam.segType = 'computed';
    % visualizeSegFamily(computedSegments,paramVisSegFam);
    % title('Computed segmentation');
    %}

    % attach audio file to segment family visualization
    %{
    % parameterMPP.featureRate = parameter.featureRate;
    % parameterMPP.fs = sideinfo.wav.fs;
    % parameterMPP.featureTimeResType = 'seconds';
    % h_fig = gcf;
    % makePlotPlayable(f_audio, h_fig, parameterMPP);
    % by left clicking on the x-axis of the figure, the playback will
    % jump to the clicked position.
    % by right clicking on the x-axis of the figure, the player will stop.
    %}

%     eval(['save' parameter.fitnessSaveFileName(1:end-1) ' S_final fitness_info parameter pathFamily']);
 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   7. Loads ground truth segmentation and compares with computed 
    %      segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % reading ground truth from txt file:

    dirAnnotation       = '../data_annotation/';
    parameter.title     = filename;
    groundTruth_struct  = structure_annotation([dirAnnotation 'pei_anno/' filename '.txt'], 1);%parseAnnotationFile([dirAnnotation parameter.title '_bar.txt']);

    % paramVisSegFam.segType = 'groundtruth';
    % visualizeSegFamily(groundTruth_struct,paramVisSegFam);
    % title('Ground truth segmentation');
    % h_fig = gcf;
    % makePlotPlayable(f_audio, h_fig, parameterMPP);

    for i=1:size(computedSegments, 2)
        Compute(computedSegments(i).start:computedSegments(i).end,1) = 1;
    end

    % show ground truth and computed result together
    figure;
    h_fig = subplot(2,1,1);
    paramVisSegFam.segType = 'groundtruth';
    visualizeSegFamily(groundTruth_struct,paramVisSegFam,h_fig);
    title('Ground truth segmentation');
    h_fig = subplot(2,1,2);
    paramVisSegFam.segType = 'computed';
    visualizeSegFamily(computedSegments,paramVisSegFam,h_fig);
    title('Computed segmentation');
    % makePlotPlayable(f_audio, h_fig, parameterMPP);

    saveas(h_fig, ['data_image/'  multiProcess '/' filename '_result_BAR.jpg']);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   8. Compute ans save Precision、Recall、F-score
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ground  = zeros(size(S_final,1), 1);
    Compute = zeros(size(S_final,1), 1);

    for i=1:size(groundTruth_struct, 2)
        Ground (groundTruth_struct(i).start:groundTruth_struct(i).end,1) = 1;   
    end
    for i=1:size(computedSegments, 2)
        Compute(computedSegments(i).start:computedSegments(i).end,1) = 1;
    end

    Evaluation.true_positive  = sum( Ground .* Compute  );
    Evaluation.false_negtive  = sum((Ground-Compute)== 1);
    Evaluation.false_positive = sum((Ground-Compute)==-1);

    Evaluation.Precision = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_positive) , 4);
    Evaluation.Recall    = round(Evaluation.true_positive/(Evaluation.true_positive+Evaluation.false_negtive)  , 4);
    Evaluation.F_score   = round(2*Evaluation.Precision *Evaluation.Recall/(Evaluation.Precision+Evaluation.Recall), 4);
    eval(['save ../evaluation/'  multiProcess '/' filename '_Coarse_result_' num2str(S_no) ' computedSegments groundTruth_struct Evaluation']);

    close  all;
    
    %% 精細分析
%     紀錄pathFamily每個對應的bar的相似度從哪個similar matrix來的

    % 找最後的pathfamily
    p = 1;
    maxValue = 0;
    for k = 1:length(DTWpathFamily)
        if mean(mean(DTWpathFamily{k,1})) > maxValue
            maxValue = mean(mean(DTWpathFamily{k,1}));
            p = k;
        end
    end
    
%     for p = 1:length(DTWpathFamily)
        clear PATHFAMILY p1 p1 ; close all;
        for j=1:length(DTWpathFamily{p,1})
            if strcmp(multiProcess, 'allPath')
                DTWpathFamily{p,1}(3,j) = S.index(DTWpathFamily{p,1}(1,j),DTWpathFamily{p,1}(2,j));

                % 比較路徑 p1 p2從哪來 right left or combine
                if DTWpathFamily{p,1}(3,j)==1
                    PATHFAMILY(p).p1(j,1) = P_right.down_PATH_Y( DTWpathFamily{p,1}(1,j), 1 );
                    PATHFAMILY(p).p2(j,1) = P_right.down_PATH_Y( DTWpathFamily{p,1}(2,j), 1 );
                elseif DTWpathFamily{p,1}(3,j)==2
                    PATHFAMILY(p).p1(j,1) = P_left.down_PATH_Y ( DTWpathFamily{p,1}(1,j), 1 );
                    PATHFAMILY(p).p2(j,1) = P_left.down_PATH_Y ( DTWpathFamily{p,1}(2,j), 1 );
                elseif DTWpathFamily{p,1}(3,j)==3
                    PATHFAMILY(p).p1(j,1) = P_right.down_PATH_Y( DTWpathFamily{p,1}(1,j), 1 );
                    PATHFAMILY(p).p2(j,1) = P_left.down_PATH_Y ( DTWpathFamily{p,1}(2,j), 1 );
                end
                
            else
                DTWpathFamily{p,1}(3,j) = S.index(DTWpathFamily{p,1}(1,j),DTWpathFamily{p,1}(2,j));
                % 比較路徑 p1 p2從哪來 right left or combine
                if DTWpathFamily{p,1}(1,end)>size(P_right.down_PATH_Y,2); DTWpathFamily{p,1}(1,end) = size(P_right.down_PATH_Y,2);  end
                if DTWpathFamily{p,1}(3,j)==1
                    p1(:,j) = P_right.down_PATH_Y( :, DTWpathFamily{p,1}(1,j) );
                    p2(:,j) = P_right.down_PATH_Y( :, DTWpathFamily{p,1}(2,j) );
                elseif DTWpathFamily{p,1}(3,j)==2
                    p1(:,j) = P_left.down_PATH_Y ( :, DTWpathFamily{p,1}(1,j) );
                    p2(:,j) = P_left.down_PATH_Y ( :, DTWpathFamily{p,1}(2,j) );
                elseif DTWpathFamily{p,1}(3,j)==3
                    p1(:,j) = P_right.down_PATH_Y( :, DTWpathFamily{p,1}(1,j) );
                    p2(:,j) = P_left.down_PATH_Y ( :, DTWpathFamily{p,1}(2,j) );
                end
                
            end
        end
        [ PATHFAMILY, SHIFTMATRIX ]        = ER_key_relation( p1, p2, DTWpathFamily{p,1}, filename, multiProcess, p, S_no );
%         if p==1
%             if strcmp(multiProcess, 'allPath')
%                 [ PATHFAMILY_v2,SHIFTMATRIX ]      = ER_key_relation_v2( PATHFAMILY(p).p1, PATHFAMILY(p).p2, DTWpathFamily{p,1}, filename, multiProcess );
%             else
%                 [ PATHFAMILY, SHIFTMATRIX ]        = ER_key_relation( p1, p2, DTWpathFamily{p,1}, filename, multiProcess);
%             end
%         end

%     end

end

end  
end    
% end
    