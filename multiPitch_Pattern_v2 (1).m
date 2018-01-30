function [ result ] = multiPitch_Pattern_v2( midiInfo, timeSig, plotPattern )

%%�e�X�C�p�`�� pattern 
%   input  : midiInfo, case, plot
%   output : �C�p�`pattern(�p�`�I��*�p�`��)

 %%%%%%%%%%%%%%%%% �Y�Ϥ��P���縹�A�]�i�H�e�� (�h���B�z)
 
    if nargin < 3, plotPattern  = 0;  end

%% testing code
%     clear all; close all; clc;
%     [ midiInfo, timeSig ] = midi_Preprocess('../midi/b_8_1', 1);    % �ݭn midiinfo 1(onset beat) 4(pitch)
%     plotPattern = 1;

    
%% initial
    Inner_point     = 100;
    resultPattern   =  [];
    patternPoint    =  50;
    addBar          =   0;
    addBeat         =   0;
    
%%
    for i = 1:size(timeSig,1)
        % set
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   ceil( (timeSig(i+1,5) - timeSig(i,5) ) /NoOfBeats);
        else
            NoOfBar     =   ceil( ( midiInfo(end,1) + midiInfo(end,2) - timeSig(i,5) ) / NoOfBeats);
        end
        
        for j = 1:NoOfBar
            
            note_X = []; note_Y = [];
            bar_start   = (addBeat+(j-1)*NoOfBeats);
            bar_end     = addBeat+j*NoOfBeats;
            notesOfBar  = intersect(find(midiInfo(:,1)>=bar_start), find(midiInfo(:,1)<bar_end));

            if ~isempty(notesOfBar) % ���p�`��������
                tab         =   tabulate(midiInfo(notesOfBar,1));
                tab(tab(:,2)==0,:) = [];
                
                % �p�`��T
                Bar(addBar+j).notePos(:,1)     = tab(:,1); % �X�Ӧ�m��note
                Bar(addBar+j).eachPosNum(:,1)  = tab(:,2); % �C�Ӧ�mnote�ƶq
                for k=1:length(tab(:,1)) % �C�Ӧ�mnote pitch number(�s��)
                    Bar(addBar+j).pitch{:,k}   = midiInfo(notesOfBar(midiInfo(notesOfBar,1) == tab(k,1)), 4);   
                end         
                Bar(addBar+j).noteNum          = length(notesOfBar); % �p�`��note�`��
                Bar(addBar+j).pathNum          = prod(tab(:,2));     % �p�`��path�`��
                
                                %{
                % ����
%                 if rest
%                     % case 1 : �}�Y���
%                     interval = 0;
%                     if notesOfBar(1)~=1
%                         interval = max(bar_start, sum(midiInfo(notesOfBar(1)-1,1:2)));
%                     end
%                     
%                     if Bar(addBar+j).notePos(1,1) - interval > ((1/2)^4)
%                         Bar(addBar+j).notePos       = [bar_start; Bar(addBar+j).notePos];
%                         Bar(addBar+j).eachPosNum    = [1; Bar(addBar+j).eachPosNum     ];
%                         Bar(addBar+j).pitch         = [{eps+1}, Bar(addBar+j).pitch    ];
%                         a = 1;
%                     end
%                     % case 2 : �������
%                     rest_on = (midiInfo(notesOfBar(end),1)+midiInfo(notesOfBar(end),2));
%                     if bar_end - rest_on > (1/2)^4
%                         Bar(addBar+j).notePos       = [Bar(addBar+j).notePos; rest_on+0.01];
%                         Bar(addBar+j).eachPosNum    = [Bar(addBar+j).eachPosNum; 1        ];
%                         Bar(addBar+j).pitch         = [Bar(addBar+j).pitch, {eps+1}       ];
%                         a = 2;
%                     end
%                     % case 3 : �����������
%                     X_on  = midiInfo(notesOfBar,1);
%                     X_off = midiInfo(notesOfBar,1) + midiInfo(notesOfBar,2);
%                     find(X_on(2:end) - X_off(1:end-1) > (1/2)^4, 1);
%                     
%                     if ~isempty(find(X_on(2:end) - X_off(1:end-1) > (1/2)^4, 1)) % , 1
%                         rest_num = length(find(X_on(2:end) - X_off(1:end-1) > (1/2)^4));
%                         Bar(addBar+j).notePos        = [ Bar(addBar+j).notePos; round(X_off(X_on(2:end) - X_off(1:end-1) > (1/2)^4),2)+0.01 ];
%                         [Bar(addBar+j).notePos, ind] = sort(Bar(addBar+j).notePos);
%                         Bar(addBar+j).eachPosNum     = [Bar(addBar+j).eachPosNum; ones(rest_num, 1) ];
%                         Bar(addBar+j).eachPosNum     = Bar(addBar+j).eachPosNum(ind);
%                         Bar(addBar+j).pitch          = [Bar(addBar+j).pitch, repmat({eps+1},1,rest_num)];
%                         Bar(addBar+j).pitch          = Bar(addBar+j).pitch(ind);
%                         a = 3;                       
%                     end     
%                 end
                %}
                % �e�Ҧ�bar path
                PATH_Y = [];

                for posNum = 1:length(Bar(addBar+j).eachPosNum) % ��p�`���Ҧ����|��Y                  
                    eachPosNum = Bar(addBar+j).eachPosNum(posNum);                    
                    if posNum == 1
                        PATH_Y = Bar(addBar+j).pitch{posNum};
                    else
                        PATH_Y = repmat(PATH_Y,eachPosNum,1);
                        PATH_Y = [PATH_Y sort(repmat(Bar(addBar+j).pitch{posNum}, size(PATH_Y,1)/eachPosNum, 1))];
                    end                   
                end

                BAR(addBar+j).pathNum = prod(tab(:,2));   % �p�`��path�`��
                BAR(addBar+j).PATH_X  = (Bar(addBar+j).notePos(:,1)-(addBeat+(j-1)*NoOfBeats)) * Inner_point;
                BAR(addBar+j).PATH_Y  = PATH_Y';

            else % �p�`���S����    %%%%%% ���żȮ��ܬ��@���u pattern
                BAR(addBar+j).pathNum = 1;
                BAR(addBar+j).PATH_X  = 0;
                BAR(addBar+j).PATH_Y  = eps+1; % 0 
            end
               
             %% ���t                         
            X = BAR(addBar+j).PATH_X;
            Y = BAR(addBar+j).PATH_Y;
            
            X(X==0) = 1;                                       % matlab index �q 1 �}�l
            if X(1)~=1; X = [1; X];  Y = [Y(1,:); Y]; end      % �ɤp�`�}�Y
            Y = [Y; Y(end,:)]; X = [X; NoOfBeats*Inner_point]; % �ɤp�`����
            

                BAR(addBar+j).in_PATH_Y = interp1(X,Y,1:NoOfBeats*Inner_point);
                if size(BAR(addBar+j).in_PATH_Y, 1) ~= NoOfBeats*Inner_point
                    BAR(addBar+j).in_PATH_Y = BAR(addBar+j).in_PATH_Y';
                end
                
            if length(BAR(addBar+j).in_PATH_Y) < patternPoint
                tmp = 1/(ceil(patternPoint/length(BAR(addBar+j).in_PATH_Y)+1));
                BAR(addBar+j).in_PATH_Y = BAR(addBar+j).in_PATH_Y(1:tmp:end);
            end
            %{
         %% ���t                         
%             X = BAR(addBar+j).PATH_X;
%             Y = BAR(addBar+j).PATH_Y;
%             
%             X(X==0) = 1;                                       % matlab index �q 1 �}�l
%             if X(1)~=1; X = [1; X];  Y = [Y(1,:); Y]; end      % �ɤp�`�}�Y
%             Y = [Y; Y(end,:)]; X = [X; NoOfBeats*Inner_point]; % �ɤp�`����
%             
%             for pathInx = 1:size(Y,2)
%                 
%                 inner_Y = []; 
% 
%                 for k=2:length(X)
%                     interval = Y(k,pathInx) - Y(k-1,pathInx);
% 
%                     if interval == 0
%                         inner_Y(X(k-1):X(k), pathInx) = Y(k,pathInx);
%                     else
%                         inner_Y(X(k-1):X(k), pathInx) = linspace(Y(k-1,pathInx), Y(k,pathInx), X(k)-X(k-1)+1);
%                     end
%                 end
% 
%                 BAR(addBar+j).in_PATH_Y(:,pathInx) = inner_Y(:,pathInx);
%              end
%}
        %% downsampling
            downSample   = size(BAR(addBar+j).in_PATH_Y,1)/patternPoint;
            BAR(addBar+j).down_PATH_Y          = BAR(addBar+j).in_PATH_Y(1:downSample:end, :);    % ĵ�i�I �]�� downSample ���O���
            
            BAR(addBar+j).down_PATH_Y(end+1,:) = BAR(addBar+j).in_PATH_Y(end, :); 

         %% �p�Gpattern���O�@�˪��� corr �|��XNAN 
                % �̫�@�ӭȴ�0.001
            for pathInx = 1:size(Y,2)       
                if all( abs(diff( BAR(addBar+j).down_PATH_Y(:,pathInx) )) < 0.0001 )
                    if sum(BAR(addBar+j).down_PATH_Y(:,pathInx))<=size(BAR(addBar+j).down_PATH_Y(:,pathInx),1)+1 % ��Ӥp�`���
                        BAR(addBar+j).down_PATH_Y(end, pathInx) = BAR(addBar+j).down_PATH_Y(end-1, pathInx) + 0.0001;
                    else % ��Ӥp�`���O�@�˪����Ϊ̥u���@�ӭ�
                        BAR(addBar+j).down_PATH_Y(end, pathInx) = BAR(addBar+j).down_PATH_Y(end-1, pathInx) - 0.0001;
                    end                
                end
            end
            
        %% normalize 
%             BAR(addBar+j).nor_PATH_Y = BAR(addBar+j).down_PATH_Y./max(BAR(addBar+j).down_PATH_Y); % �C��path�ӧOnormalize
            BAR(addBar+j).nor_PATH_Y = BAR(addBar+j).down_PATH_Y/max(max(BAR(addBar+j).down_PATH_Y)); % �@�_normalize
%             NAN �G ��0�����Y
            result.barPath{addBar+j,1}     =  BAR(addBar+j).nor_PATH_Y;
            result.maxPitch(addBar+j,1)    =  max(max(BAR(addBar+j).down_PATH_Y));
            result.down_PATH_Y{addBar+j,1} =  BAR(addBar+j).down_PATH_Y;
            
        %% plot
            if plotPattern
                index = mod(j,16);
                if index==1; figure;    end
                if index==0; index=16;  end
                plot_x = linspace(0, NoOfBeats, patternPoint+1);
                note_Y = midiInfo(notesOfBar,4)./max(max(BAR(addBar+j).down_PATH_Y));
                note_X = midiInfo(notesOfBar,1)-(addBeat+(j-1)*NoOfBeats);
                subplot(4,4,index); 
                plot(plot_x, BAR(addBar+j).nor_PATH_Y); hold all; plot(note_X, note_Y, '.','MarkerSize',25);
                title(['bar ' num2str(addBar+j)]);  axis([0,NoOfBeats+1,0,1.1]); %NoOfBeats
            end 
        end

        addBar  = addBar + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
     
    
end
    