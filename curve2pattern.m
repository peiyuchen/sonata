function [ result ] = curve2pattern( midiINFO, curve, timeSig, plotPattern )

%%�e�X�C�p�`�� pattern 
%   input  : midiINFO, case, plot
%   output : �C�p�`pattern(�p�`�I��*�p�`��)

 %%%%%%%%%%%%%%%%% �Y�Ϥ��P���縹�A�]�i�H�e��
 
    if nargin < 4, plotPattern  = 0;  end

    %% initial
    innerPoint      = 100;
    resultPattern   = [];
    patternPoint    = 50;
    addBar          = 0;
    addBeat         = 0;
    
    %% pitch curve

    for i = 1:size(timeSig,1)
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   max(ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) / NoOfBeats), 1);
        else
            NoOfBar     =   max(ceil( ( length(curve)/innerPoint - timeSig(i,5) ) / NoOfBeats), 1);
        end
        normal_MAX      =   [];
        
        for j = 1:NoOfBar
            X = []; Y = []; note_X = []; note_Y = [];
            barStart   = (addBeat+(j-1)*NoOfBeats);
            barEnd     = addBeat+j*NoOfBeats;
            notesOfBar = intersect(find(midiINFO(:,1)>=barStart), find(midiINFO(:,1)<barEnd));
            currentBar = addBar + j;
            
            if ~isempty(notesOfBar) % ���p�`��������
                X = midiINFO(notesOfBar,1);
                note_Y = curve(round(X * innerPoint));
                note_X = X - (j-1)*NoOfBeats - addBeat;
            end
            
            X = [barStart+0.01; X; barEnd-0.01];
            X = round(X * innerPoint);
            X(X==0) = 1; X = unique(X);
            Y       = curve(X);

            % �e��O���ŴN�����̫�e�����A���n�O����
            if Y(1) == 0;   Y(1)=Y(2);         end  
            if Y(end) == 0; Y(end) = Y(end-1); end
      
            X = X - (j-1)*NoOfBeats*innerPoint-addBeat*innerPoint;  % �u�ݦ��p�`����m
            

            %% ���t        
            inner_Y = [];
            
            xi = linspace(1, NoOfBeats*innerPoint-1, patternPoint); % ����1�|inner_Y�|NaN
            inner_Y = interp1(X, Y, xi, 'pchip');                   % ��k�Glinear, pchip  

            barPattern(:,j) = inner_Y;
            
        %% �p�Gpattern���O�@�˪��� corr �|��XNAN 
        % �̫�@�ӭȴ�0.001
            smallValue = 10^(-3);
            if all( abs(diff(barPattern(:,j))) < smallValue )
                if sum(barPattern(:,j)) <= size(barPattern,1)+1 % ��Ӥp�`���
                    barPattern(end,j) = max(0, barPattern(end-1,j) + smallValue); 
                else % ��Ӥp�`���O�@�˪����Ϊ̥u���@�ӭ�
                    barPattern(end,j) = max(0, barPattern(end-1,j) - smallValue); 
                end
            end

            
        %% normalize 
            normal_MAX(:,j) = barPattern(:,j)/max(barPattern(:,j));

        %% plot
%             if plotPattern
%                 index = mod(j,16);
%                 if index==1; figure;    end
%                 if index==0; index=16;  end
%                 plot_x = linspace(0, NoOfBeats, patternPoint);
%                 note_Y = note_Y / max(barPattern(:,j));
%                 subplot(4,4,index); 
%                 plot(plot_x, normal_MAX(:,j)); hold all; plot(note_X, note_Y, '.','MarkerSize',25);
%                 title(['bar ' num2str(addBar+j)]);  axis([0,NoOfBeats+1,0,1.1]); %NoOfBeats
%             end   
if plotPattern
            if currentBar==196 || currentBar==376
                currentBar
                if currentBar==196; color = 'red';
                else
                    color = 'blue';
                end
                plot_x = linspace(0, NoOfBeats*2, patternPoint);
                
                plot(plot_x, barPattern(:,j),'color',color,'LineWidth',2); hold all;
                axis([-inf, inf, 50, 75]);
                xlabel('Beat'); ylabel('Pitch number');
                if currentBar==376; legend('196 bar','376 bar'); end
%                 set(gca,'FontName','Times New Roman','FontSize',50) % �]�m���жb��צr��W�١A�j�p
%                 xt = get(gca, 'XTick');
%                 set(gca, 'FontSize', 40, 'FontName','Times New Roman')
            end
end
            result.barPath(:,addBar+j)     = normal_MAX(:,j);
            result.down_PATH_Y(:,addBar+j) = barPattern(:,j);
        end
        
        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;
        resultPattern = [ resultPattern, normal_MAX ];
        
    end
end
    
