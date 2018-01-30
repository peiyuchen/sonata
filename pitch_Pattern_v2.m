% function [ resultPattern ] = pitch_Pattern_v2( midiInfo, plotPattern )

%%�e�X�C�p�`�� pattern 
%   input  : midiInfo, case, plot
%   output : �C�p�`pattern(�p�`�I��*�p�`��)

 %%%%%%%%%%%%%%%%% �Y�Ϥ��P���縹�A�]�i�H�e��
 
%     if nargin < 2, plotPattern  = 0;  end

%% testing code
%     clear all; close all; clc;

    [midiInfo,timeSig]    = midi_Preprocess('../midi/b_1_1', 2);   % �ݭn midiinfo 1(onset beat) 4(pitch)
    plotPattern = 1;
    
    
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
    curve = [];
    
    % pitch curve
        curve = [(eps+1)*ones(round(midiInfo(end,1)+midiInfo(end,2))*Inner_point, 1)];

    for j = 1:size(midiInfo,1)
        Start = round(midiInfo(j,1)*Inner_point); if Start==0; Start=1; end;
        curve( Start:round((midiInfo(j,1) + midiInfo(j,2))*Inner_point),1 ) = midiInfo(j,4);
    end

    
    for i=1:size(timeChan,1)
        % set
        NoOfBeats       =   timeChan(i,3)*(4/(timeChan(i,4)))
        NoOfBar         =   floor((midiInfo(timeChan(i,2),1)-midiInfo(timeChan(i,1),1)+midiInfo(timeChan(i,2),2))/NoOfBeats)+1 
        normal_MAX      =   [];
        for j = 1:NoOfBar
            X = []; Y = []; 
            bar_start   = (addBeat+(j-1)*NoOfBeats);
            bar_end     = addBeat+j*NoOfBeats;
            notesOfBar  = intersect(find(midiInfo(:,1)>=bar_start), find(midiInfo(:,1)<bar_end));

            if ~isempty(notesOfBar) % ���p�`��������
%                 X = round(midiInfo(notesOfBar,1)*Inner_point);
                X = midiInfo(notesOfBar,1);
          
            else % �p�`���S����    %%%%%% ���żȮ��ܬ��@���upattern
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

            %% ���t
            
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

%             patternPoint = 50;%NoOfBeats*10; % �q���� sampling �I
            downSample   = size(inner_Y,1)/patternPoint;

            if length(inner_Y)>1
                barPattern(1:patternPoint, j)   = inner_Y(1:downSample:end, 1);    % ĵ�i�I �]�� downSample ���O���
                barPattern(patternPoint+1, j)   = inner_Y(end, 1); 
            else
                barPattern(1:patternPoint+1, j) = inner_Y;
            end
            
        %% �p�Gpattern���O�@�˪��� corr �|��XNAN 
        % �̫�@�ӭȴ�0.001
            if all(abs(diff(barPattern(:,j)))<0.0001)
                barPattern(end,j) = max(0,barPattern(end-1,j)-0.0001);
            end
        %% normalize 
            normal_MAX(:,j)     = barPattern(:,j)/max(barPattern(:,j));
        
        %% plot
            if plotPattern
                index = mod(j,16);
                if index==1; figure;    end
                if index==0; index=16;  end
                plot_x = linspace(0, NoOfBeats, patternPoint+1);

                subplot(4,4,index); 
                plot(plot_x, normal_MAX(:,j));
                title(['bar ' num2str(addBar+j)]);  %axis([0,NoOfBeats+1,0.8,1.01]); %NoOfBeats
            end    
            
        end
        
        addBar  = addBar + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;
        resultPattern = [resultPattern, normal_MAX ];
        
    end
% end
    
