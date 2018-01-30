% function RESULT = structure_annotation(filename, type)
filename = '../data_annotation/paper_anno/b_1_1.txt'; type = 1;
% 把paper標記的變成我要的格式

% http://webcache.googleusercontent.com/search?q=cache:WCBUkOB0oeAJ:www.ecaa.ntu.edu.tw/weifang/matlab%2520tutor/18-%25E6%25AA%2594%25E6%25A1%2588%25E8%25AE%2580%25E5%25AF%25AB.ppt+&cd=9&hl=zh-TW&ct=clnk&gl=tw&client=safari

%     dir = '../data_annotation/paper_anno/'; filename = [dir 'b_1_1.meta'];

%     dir = '../data_annotation/pei_anno/'; filename = [dir 'b_16_1.txt'];
        
%     type = 1;

    fid = fopen(filename, 'r');

    i = 1; j = 1; E_no = 1; R_no = 1; row_no = 1; ROW_NO = 1; coarse_end_index = [];

    while feof(fid) == 0 % feof 測試檔案指標是否已到達結束位置

          line{i,1} = fgetl(fid);

%           disp(line);

          if length(line{i,1}) > 9

              if strcmp(line{i,1}(1:9), 'structure')	% 比較兩字串 str1 和 str2

                  data_anno(j,1:5) = {[],[],[],[],[]};

                  tmp  = regexp(line{i,1}, '\s+', 'split');

                  tmpp = regexp(tmp{4}, '"+', 'split');

                  tmppp = regexp(tmpp{2}, ':', 'split');

                  tmpppp = [tmp(1:3) tmppp];

                  data_anno(j,1:length(tmpppp)) = tmpppp;
                  
                  data_anno(j,2) = { ceil( str2num( data_anno{j,2} ) ) };
                  
                  if data_anno{j,2} == 0; data_anno{j,2} = 1; end
                  
                  if floor( str2num( data_anno{j,3} ) + 0.001 ) > floor( str2num( data_anno{j,3} ) )
                      
                      data_anno(j,3) = { floor( str2num( data_anno{j,3} ) ) };
                  
                  else
                      
                      data_anno(j,3) = { floor( str2num( data_anno{j,3} ) ) + 1 };
                      
                  end

                  % Exposition

                  if strncmp(data_anno{j,4}, 'Exposition', 10) && strncmp(data_anno{j,5}, 'FirstGroup', 10)     

%                       result(row_no, 1:2) = data_anno(j, 2:3);

%                       result(row_no, 5)   = {['E' num2str(E_no)]};
                      
                      GT_fine(row_no).M1_start = data_anno{j, 2};
                      
                      GT_fine(row_no).M1_end = data_anno{j, 3};
                      
                      GT_coarse(ROW_NO).start = data_anno{j, 2};
                      
                      GT_coarse(ROW_NO).end = [];
                      
                      GT_coarse(ROW_NO).label = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];
                      
                      ROW_NO = ROW_NO + 1;

                  end

                  if strncmp(data_anno{j,4}, 'Exposition', 10) && strncmp(data_anno{j,5}, 'SecondGroup', 11)

%                       result(E_no, 3:4) = data_anno(j, 2:3);
                      
                      GT_fine(row_no).M2_start = data_anno{j, 2};
                      
                      GT_fine(row_no).M2_end = data_anno{j, 3};
                      
                      GT_fine(row_no).label = data_anno{j, 4};

                      E_no = E_no + 1;    

                      row_no = row_no + 1;   

                  end
                  
                  % Development
                  if ~strncmp(data_anno{j,4}, 'Exposition', 10) && ~strncmp(data_anno{j,4}, 'Recapitulation', 14)
                      
                      GT_coarse(ROW_NO).start = data_anno{j, 2};
                      
                      GT_coarse(ROW_NO).end = [];
                      
                      GT_coarse(ROW_NO).label = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];
                      
                      ROW_NO = ROW_NO + 1;                  
                  
                  end
                  
                  % Recapitulation

                  if strncmp(data_anno{j,4}, 'Recapitulation', 14) && strncmp(data_anno{j,5}, 'FirstGroup', 10)

%                       result(row_no, 1:2) = data_anno(j, 2:3);

%                       result(row_no, 5)   = {['R' num2str(R_no)]};
                      
                      GT_fine(row_no).M1_start = data_anno{j, 2};
                      
                      GT_fine(row_no).M1_end = data_anno{j, 3};
                      
                      
                      GT_coarse(ROW_NO).start = data_anno{j, 2};
                      
                      GT_coarse(ROW_NO).end = [];
                      
                      GT_coarse(ROW_NO).label = data_anno{j, 4};
                      
                      coarse_end_index = [coarse_end_index j];
                      
                      ROW_NO = ROW_NO + 1;

                  end

                  if strncmp(data_anno{j,4}, 'Recapitulation', 14) && strncmp(data_anno{j,5}, 'SecondGroup', 11)

%                       result(row_no, 3:4) = data_anno(j, 2:3);
                      
                      GT_fine(row_no).M2_start = data_anno{j, 2};
                      
                      GT_fine(row_no).M2_end = data_anno{j, 3};
                      
                      GT_fine(row_no).label = data_anno{j, 4};

                      R_no = R_no + 1;   

                      row_no = row_no + 1;
                      

                  end


                  j = j + 1;

                  tmp = []; tmpp = []; tmppp = []; tmpppp = [];

              end

          end

          i = i + 1;

    end
    
    
    coarse_end_index(1) = [];
    
    k = 1;
    
    for i = 1 : length(coarse_end_index)
    
        GT_coarse(i).end = data_anno{coarse_end_index(i)-1, 3};
        
        if ~strncmp(GT_coarse(i).label, 'Development', 11)
            
            GT_coarse_ER(k).start = GT_coarse(i).start;
            
            GT_coarse_ER(k).end = GT_coarse(i).end;
        
            GT_coarse_ER(k).label = GT_coarse(i).label;
            
            k = k + 1;
            
        end
   
    end
    
    GT_coarse(length(coarse_end_index)+1).end = data_anno{size(data_anno, 1), 3};
    
    if strncmp(GT_coarse(length(coarse_end_index)+1).label, 'Recapitulation', 11)
        
            GT_coarse(length(coarse_end_index)+1).end = data_anno{size(data_anno, 1), 3};
    
            GT_coarse_ER(k).start = GT_coarse(end).start;
            
            GT_coarse_ER(k).end = GT_coarse(end).end;
        
            GT_coarse_ER(k).label = GT_coarse(end).label;
    
    end
    
    if type == 1
        
        RESULT = GT_coarse_ER;
    
    elseif type == 2
    
        RESULT = GT_fine;
    
    end
      

end










