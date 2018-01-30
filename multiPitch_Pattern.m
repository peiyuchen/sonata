function [ result ] = multiPitch_Pattern( midiInfo, timeSig, plotPattern )

%%畫出每小節的 pattern 
%   input  : midiInfo, case, plot
%   output : 每小節pattern(小節點數*小節數)

 %%%%%%%%%%%%%%%%% 即使不同的拍號，也可以畫圖 (多音處理)
 
    if nargin < 3, plotPattern  = 0;  end

%% testing code
%     clear all; close all; clc;
%     midiInfo    = midi_Preprocess('../midi/b_1_1', 1);    % 需要 midiinfo 1(onset beat) 4(pitch)
%     plotPattern = 0;

    
%% initial
    Inner_point     =   100;
%     C4              =   60;
    resultPattern   = [];
    patternPoint    = 50;
    addBar          = 0;
    addBeat         = 0;
    
%% time Sig 

    timeChan(:,2)       = find(diff(midiInfo(:,9))~=0); timeChan(end+1,2) = size(midiInfo,1);
    timeChan(2:end,1)   = timeChan(1:end-1,2)+1;        timeChan(1,1) = 1;
    timeChan(:,3)       = midiInfo(timeChan(:,1), 9);
    timeChan(:,4)       = midiInfo(timeChan(:,1),10);
    
%% pitch curve

    for i=1:size(timeSig,1)
        % set
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   ceil( (timeSig(i+1,5) - timeSig(i,5) ) /NoOfBeats);
        else
            NoOfBar     =   ceil( ( midiInfo(end,1) + midiInfo(end,2) - timeSig(i,5) ) / NoOfBeats);
        end
        
        for j = 1:NoOfBar
            
            notesOfBar  =   intersect(find(midiInfo(:,1)>=(addBeat+(j-1)*NoOfBeats)), find(midiInfo(:,1)<addBeat+j*NoOfBeats));
            
            if ~isempty(notesOfBar) % 此小節內有音符
                tab         =   tabulate(midiInfo(notesOfBar,1));
                tab(tab(:,2)==0,:) = [];
                
                % 小節資訊
                Bar(addBar+j).notePos(:,1)     = tab(:,1); 
                Bar(addBar+j).noteCount(:,1)   = tab(:,2);
                for k=1:length(tab(:,1))
                    Bar(addBar+j).pitch{:,k}   = midiInfo(notesOfBar(midiInfo(notesOfBar,1) == tab(k,1)), 4);    
                end                
                Bar(addBar+j).noteNum          = length(notesOfBar);
                Bar(addBar+j).pathNum          = prod(tab(:,2));
                
                % 畫所有bar path
                path_X = []; PATH_Y = [];
                
                for posNum = 1:length(Bar(addBar+j).noteCount)  % 一小節幾個位置
                    
                    noteCount = Bar(addBar+j).noteCount(posNum);
                    
                    if posNum == 1
                        PATH_Y = Bar(addBar+j).pitch{posNum};
                    else
                        PATH_Y = repmat(PATH_Y,noteCount,1);
                        PATH_Y = [PATH_Y sort(repmat(Bar(addBar+j).pitch{posNum}, size(PATH_Y,1)/noteCount, 1))];
                    end
                    
                end
                
                BAR(addBar+j).pathNum = prod(tab(:,2));
                BAR(addBar+j).PATH_X  = (Bar(addBar+j).notePos(:,1)-(addBeat+(j-1)*NoOfBeats)) * Inner_point;
                BAR(addBar+j).PATH_Y  = PATH_Y';
           
            else % 小節內沒音符    %%%%%% 休止符暫時變為一條線 pattern
                BAR(addBar+j).pathNum = 1;
                BAR(addBar+j).PATH_X  = 0;
                BAR(addBar+j).PATH_Y  = eps+1; % 0 
            end
                        
        %% 內差                 
            
            X = BAR(addBar+j).PATH_X;
            Y = BAR(addBar+j).PATH_Y;
            X(X==0) = 1;                                       % matlab index 從 1 開始
            if X(1)~=1; X = [1; X];  Y = [Y(1,:); Y]; end      % 補小節開頭
            Y = [Y; Y(end,:)]; X = [X; NoOfBeats*Inner_point]; % 補小節結尾
            
            for pathInx = 1:size(Y,2)
                
                inner_Y = []; 

                for k=2:length(X)
                    interval = Y(k,pathInx) - Y(k-1,pathInx);

                    if interval == 0
                        inner_Y(X(k-1):X(k), pathInx) = Y(k);
                    else
                        inner_Y(X(k-1):X(k), pathInx) = linspace(Y(k-1,pathInx), Y(k,pathInx), X(k)-X(k-1)+1);
                    end
                end

                BAR(addBar+j).in_PATH_Y(:,pathInx) = inner_Y(:,pathInx);
            end
            
        %% downsampling

            downSample   = size(BAR(addBar+j).in_PATH_Y,1)/patternPoint;

            BAR(addBar+j).down_PATH_Y          = BAR(addBar+j).in_PATH_Y(1:downSample:end, :);    % 警告！ 因為 downSample 不是整數
            BAR(addBar+j).down_PATH_Y(end+1,:) = BAR(addBar+j).in_PATH_Y(end, :); 

         %% 如果pattern都是一樣的值 corr 會算出NAN 
                % 最後一個值減0.001
            for pathInx = 1:size(Y,2)       
                if all( abs(diff( BAR(addBar+j).down_PATH_Y(:,pathInx) )) < 0.0001 )
                    BAR(addBar+j).down_PATH_Y(end, pathInx) = BAR(addBar+j).down_PATH_Y(end-1, pathInx) - 0.0001;
                end
            end
            
        %% normalize 
            BAR(addBar+j).nor_PATH_Y = BAR(addBar+j).down_PATH_Y./max(BAR(addBar+j).down_PATH_Y);
%             NAN ： 跟0有關係
            result.barPath{addBar+j,1} = BAR(addBar+j).nor_PATH_Y;
        end

        addBar  = addBar + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
    
%     rowCol=1;
%     Pattern1 = result;Pattern2 = result;
% %     h = waitbar(0,'Please wait...');
% %         figure;  
%     for i=1:1%length(Pattern1.barPath)
%         barPath1 = Pattern1.barPath{i,:};figure;  
%          subplot(2,1,1); plot(barPath1);
%         for j=1:3%length(Pattern2.barPath)
%             barPath2 = Pattern2.barPath{j,:};
%                      subplot(2,1,2); plot(barPath2);
%             tmp_S_matrix = similar_matrix(barPath1, barPath2, 'correlation');
%             
%             max1 = max(tmp_S_matrix,[],1);
%             max2 = max(tmp_S_matrix,[],2);
%             if rowCol == 1
%                 if length(max1) > length(max2)
%                     MAX = sum(max1)/length(max1);
%                 else
%                     MAX = sum(max2)/length(max2);
%                 end
%             else
%                 if length(max1) < length(max2)
%                     MAX = sum(max1)/length(max1);
%                 else
%                     MAX = sum(max2)/length(max2);
%                 end
%             end
%             
%             S{i,j}        = tmp_S_matrix;
%             S_matrix(i,j) = MAX;
%             
%         end
% %         waitbar(i/length(Pattern1.barPath))
%     end
% %     close(h) 
%     
%     %% plot similarity matrix
%     S_plot=1;
%     if S_plot
%         figure, imagesc(S_matrix); colorbar; title(['similar']);
%     end  
    
end
    

%     k = 6
%     figure;subplot(2,1,1); plot(BAR(k).in_PATH_Y);subplot(2,1,2); plot(BAR(k).down_PATH_Y);
% figure; plot(result.barPath{k,1});