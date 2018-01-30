function [ resultPattern ] = pitch_pattern_v3( midiInfo, timeSig, plotPattern )

%%畫出每小節的 pattern 
%   input  : midiInfo, case, plot
%   output : 每小節pattern(小節點數*小節數)

 %%%%%%%%%%%%%%%%% 即使不同的拍號，也可以畫圖
 
    if nargin < 3, plotPattern  = 0;  end

%% testing code
%     clear all; close all; clc;
% 
%     [midiInfo,timeSig]    = midi_Preprocess('../midi/b_1_1', 2);   % 需要 midiinfo 1(onset beat) 4(pitch)
%     plotPattern = 1;
    
    
%% initial
    Inner_point     =   100;
%     C4              =   60;
    resultPattern   = [];
    patternPoint    = 50;
    addBar          = 0;
    addBeat         = 0;
    
%% pitch curve

    % pitch curve
%         curve = [(eps+1)*ones(round(midiInfo(end,1)+midiInfo(end,2))*Inner_point, 1)];
    curve = [(eps+1)*ones((ceil((midiInfo(end,1)+midiInfo(end,2)-timeSig(end,5))/(timeSig(end,1)*(4/(2^timeSig(end,2)))))* ... 
                                                        (timeSig(end,1)*(4/(2^timeSig(end,2))))+timeSig(end,1))*Inner_point, 1)];


    for j = 1:size(midiInfo,1)
        Start = round(midiInfo(j,1)*Inner_point); if Start==0; Start=1; end;
        curve( Start:round((midiInfo(j,1) + midiInfo(j,2))*Inner_point),1 ) = midiInfo(j,4);
    end

    for i=1:size(timeSig,1)
        % set
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) /NoOfBeats);
        else
            NoOfBar     =   ceil( ( midiInfo(end,1) + midiInfo(end,2) - timeSig(i,5) ) / NoOfBeats);
        end
        normal_MAX      =   [];

        for j = 1:NoOfBar
            X = []; Y = []; note_X = []; note_Y = [];
            bar_start   = (addBeat+(j-1)*NoOfBeats);
            bar_end     = addBeat+j*NoOfBeats;
            notesOfBar  = intersect(find(midiInfo(:,1)>=bar_start), find(midiInfo(:,1)<bar_end));

            if ~isempty(notesOfBar) % 此小節內有音符
%                 X = round(midiInfo(notesOfBar,1)*Inner_point);
                X = midiInfo(notesOfBar,1);
                
                % 休止符
                    % case 1 : 開頭休止
                    if X(1) - bar_start > ((1/2)^4)
                        X = [bar_start; X];
                        
                        a=1;
                    end
                    % case 2 : 結束休止
                    rest_on = (X(end)+midiInfo(notesOfBar(end),2));
                    if bar_end - rest_on > (1/2)^4
                        X = [X; rest_on+0.01];
                        a=2;
                    end
                    % case 3 : 音的中間休止
                    X_on  = midiInfo(notesOfBar,1);
                    X_off = midiInfo(notesOfBar,1) + midiInfo(notesOfBar,2);
                    
                    if ~isempty(find(X_on(2:end) - X_off(1:end-1) > (1/2)^4, 1))
                        X;
                        X = [X; round(X_off(X_on(2:end) - X_off(1:end-1) > (1/2)^4),2)+0.01];
                        X = sort(X) ;
                        a = 3;
                        
                    end      
                    note_Y = curve(round(X * Inner_point));
                    note_X = X - (j-1)*NoOfBeats - addBeat;
%                     end
            else % 小節內沒音符    %%%%%% 休止符暫時變為一條線pattern
%                 X = [addBeat+(j-1)*NoOfBeats+1; addBeat+j*NoOfBeats-1]*Inner_point;
                X = [bar_start+1; bar_end-1];
            end
            
            X = round(X * Inner_point);
            X(X==0) = 1; 
            Y = curve(X);
            X = X - (j-1)*NoOfBeats*Inner_point-addBeat*Inner_point;
            Y = [Y; Y(end)]; X = [X; NoOfBeats*Inner_point];

            X(X==0) = 1; 
            if X(1)~=1; X = [1; X];  Y = [Y(1); Y]; end

            %% 內差
            
            inner_Y = [];
            if length(X)>1
                for k=2:length(X)
                    interval = Y(k)-Y(k-1);

                    if interval == 0
                        inner_Y(X(k-1):X(k), 1) = Y(k);
                    else
                        inner_Y(X(k-1):X(k), 1) = linspace(Y(k-1), Y(k), X(k)-X(k-1)+1);
                    end
                end
            else
                inner_Y = Y;
            end
            innter_Y =  curve(NoOfBeats*Inner_point*(i-1)+1:NoOfBeats*Inner_point*i);
            
        %% downsampling

%             patternPoint = 50;%NoOfBeats*10; % 訂死的 sampling 點
            downSample   = size(inner_Y,1)/patternPoint;

            if length(inner_Y)>1
                barPattern(1:patternPoint, j)   = inner_Y(1:downSample:end, 1);    % 警告！ 因為 downSample 不是整數
                barPattern(patternPoint+1, j)   = inner_Y(end, 1); 
            else
                barPattern(1:patternPoint+1, j) = inner_Y;
            end
            
        %% 如果pattern都是一樣的值 corr 會算出NAN 
        % 最後一個值減0.001
            if all(abs(diff(barPattern(:,j)))<0.0001)
                if sum(barPattern(:,j))<=size(barPattern,1)+1 % 整個小節休止
                    barPattern(end,j) = max(0,barPattern(end-1,j)+0.0001); 
                else % 整個小節都是一樣的音或者只有一個音
                    barPattern(end,j) = max(0,barPattern(end-1,j)-0.0001); 
                end
            end
        %% normalize 
            normal_MAX(:,j)     = barPattern(:,j)/max(barPattern(:,j));
        
        %% plot
            if plotPattern
                index = mod(j,16);
                if index==1; figure;    end
                if index==0; index=16;  end
                plot_x = linspace(0, NoOfBeats, patternPoint+1);
                note_Y = note_Y / max(barPattern(:,j));
                subplot(4,4,index); 
                plot(plot_x, normal_MAX(:,j)); hold all; plot(note_X, note_Y, '.','MarkerSize',25);
                title(['bar ' num2str(addBar+j)]);  axis([0,NoOfBeats+1,0,1.1]); %NoOfBeats
            end    
            
        end
        
        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;
        resultPattern = [resultPattern, normal_MAX ];
%         休止符根只有一個音差別！
    end
end
    
