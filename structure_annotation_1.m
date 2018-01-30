function RESULT = structure_annotation(filename, type)
% clear all; close; clc; filename = 'b_4_1'; type=1;

% 讀paper標記檔案
% input : 檔名, 選擇輸出類型 1.粗略結構 2.精細結構
addpath('toolbox/midi_lib/midi_lib');

[midiINFO, timeSig] = midi_Preprocess(filename);
annoFilePath        = ['../data_annotation/paper_anno/' filename '.meta'];
barAnno             = barStartEnd( midiINFO, timeSig );
fid                 = fopen(annoFilePath, 'r');

i = 1; j = 1; E_no = 1; R_no = 1; row_no = 1; ROW_NO = 1; coarse_end_index = [];

while feof(fid) == 0    % feof 測試檔案指標是否已到達結束位置

      line{i,1} = fgetl(fid);

%           disp(line);
      if length(line{i,1}) > 9
          if strcmp(line{i,1}(1:13), 'timeSignature') 
             sep        = regexp(line{i,1}, '\s+', 'split'); 
             timeSep    = [ str2num(sep{2}(1)), str2num(sep{2}(3))];
          end

          if strcmp(line{i,1}(1:9), 'structure')	% 比較兩字串 str1 和 str2

              data_anno(j,1:5) = {[],[],[],[],[]};
              
              sep            = regexp(line{i,1}, '\s+', 'split');    % 一行分成4個
              tmpSep         = regexp(   sep{4},  '"+', 'split');    % 把後面註解分開 "Exposition:FirstGroup"
              annoSep        = regexp(tmpSep{2},   ':', 'split');    % 把後面註解分開 "Exposition:FirstGroup"
%               finalSep       = [sep(1:3) annoSep];                 
              
              onset     = str2num(sep{2});
              offset    = str2num(sep{3});
              
              if j==1; onsetFirst = onset; end
              
              align1 = ( onset  - ( onsetFirst - midiINFO(1,1) ) ) * (timeSep(1)/(timeSep(2)/4));   % 對齊
              align2 = ( offset - ( onsetFirst - midiINFO(1,1) ) ) * (timeSep(1)/(timeSep(2)/4));   % 對齊
              finalOnset    = intersect(find(align1 >= barAnno(:,1)), find(align1 < barAnno(:,2)));
              finalOffset   = intersect(find(align2 >= barAnno(:,1)), find(align2 < barAnno(:,2)));
              
              if isempty(finalOffset); finalOffset = size(barAnno,1); end
              
              data_anno(j,1:3) = {sep{1}, finalOnset, finalOffset}
              data_anno(j,4:3+length(annoSep)) = annoSep;
              %%
              
               % Exposition

                  if strncmp(data_anno{j,4}, 'Exposition', 10) && strncmp(data_anno{j,5}, 'FirstGroup', 10)     

%                       result(row_no, 1:2) = data_anno(j, 2:3);
%                       result(row_no, 5)   = {['E' num2str(E_no)]};
                      
                      GT_fine(row_no).M1_start  = data_anno{j, 2};
                      GT_fine(row_no).M1_end    = data_anno{j, 3};
                      
                      GT_coarse(ROW_NO).start   = data_anno{j, 2}; 
                      GT_coarse(ROW_NO).end     = [];
                      GT_coarse(ROW_NO).label   = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];   
                      ROW_NO = ROW_NO + 1;

                  end

                  if strncmp(data_anno{j,4}, 'Exposition', 10) && strncmp(data_anno{j,5}, 'SecondGroup', 11)

%                       result(E_no, 3:4) = data_anno(j, 2:3);
                     
                      GT_fine(row_no).M2_start  = data_anno{j, 2};                     
                      GT_fine(row_no).M2_end    = data_anno{j, 3};                    
                      GT_fine(row_no).label     = data_anno{j, 4};

                      E_no   = E_no + 1;    
                      row_no = row_no + 1;   

                  end
                  
                  % Development
                  if ~strncmp(data_anno{j,4}, 'Exposition', 10) && ~strncmp(data_anno{j,4}, 'Recapitulation', 14)
                      
                      GT_coarse(ROW_NO).start   = data_anno{j, 2};                      
                      GT_coarse(ROW_NO).end     = [];
                      GT_coarse(ROW_NO).label   = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];                      
                      ROW_NO = ROW_NO + 1;                  
                  
                  end
                  
                  % Recapitulation

                  if strncmp(data_anno{j,4}, 'Recapitulation', 14) && strncmp(data_anno{j,5}, 'FirstGroup', 10)

%                       result(row_no, 1:2) = data_anno(j, 2:3);
%                       result(row_no, 5)   = {['R' num2str(R_no)]};
                      
                      GT_fine(row_no).M1_start  = data_anno{j, 2};
                      GT_fine(row_no).M1_end    = data_anno{j, 3};
                      
                      GT_coarse(ROW_NO).start   = data_anno{j, 2};
                      GT_coarse(ROW_NO).end     = [];
                      GT_coarse(ROW_NO).label   = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];
                      ROW_NO = ROW_NO + 1;

                  end

                  if strncmp(data_anno{j,4}, 'Recapitulation', 14) && strncmp(data_anno{j,5}, 'SecondGroup', 11)

%                       result(row_no, 3:4) = data_anno(j, 2:3);
                      
                      GT_fine(row_no).M2_start  = data_anno{j, 2};
                      GT_fine(row_no).M2_end    = data_anno{j, 3};
                      GT_fine(row_no).label     = data_anno{j, 4};

                      R_no   = R_no + 1;   
                      row_no = row_no + 1;

                  end           
              j = j + 1;
          end
      end
       i = i + 1;
end

coarse_end_index(1) = [];
    
    k = 1;
    
    for i = 1 : length(coarse_end_index)
    
        GT_coarse(i).end = data_anno{coarse_end_index(i)-1, 3};
        
        if ~strncmp(GT_coarse(i).label, 'Development', 11)
            
            GT_coarse_ER(k).start   = GT_coarse(i).start;
            GT_coarse_ER(k).end     = GT_coarse(i).end;
            GT_coarse_ER(k).label   = GT_coarse(i).label;
            
            k = k + 1;
            
        end
    end
    
    GT_coarse(length(coarse_end_index)+1).end = data_anno{size(data_anno, 1), 3};
    
    if strncmp(GT_coarse(length(coarse_end_index)+1).label, 'Recapitulation', 11)
        
            GT_coarse(length(coarse_end_index)+1).end = data_anno{size(data_anno, 1), 3};
            GT_coarse_ER(k).start   = GT_coarse(end).start;
            GT_coarse_ER(k).end     = GT_coarse(end).end;
            GT_coarse_ER(k).label   = GT_coarse(end).label;
    end
    
    if type == 1     
        RESULT = GT_coarse_ER;
    elseif type == 2
        RESULT = GT_fine;
    end
end
%%


function [ barAnno ] = barStartEnd( midiINFO, timeSig )
% 紀錄每個小節onset offset (beat)
    addBeat = 0; addBar = 0; midiSize = size(midiINFO,2);
    
    for i = 1:size(timeSig,1)
        NoOfBeats       =   timeSig(i,1)*(4/(2^timeSig(i,2)));
        if i ~= size(timeSig,1)
            NoOfBar     =   max(ceil( ((timeSig(i+1,5)-1) - timeSig(i,5) ) / NoOfBeats), 1);
        else
            NoOfBar     =   max(ceil( ( sum(midiINFO(end,1:2)) - timeSig(i,5) ) / NoOfBeats), 1);
        end

        for j = 1:NoOfBar
            barStart                = (addBeat+(j-1)*NoOfBeats);
            barEnd                  = addBeat+j*NoOfBeats;       
            currentBar              = addBar + j;
            barAnno(currentBar,:)   = [barStart barEnd];
        end

        addBar  = addBar  + NoOfBar;
        addBeat = addBeat + NoOfBeats*NoOfBar;

    end
end
