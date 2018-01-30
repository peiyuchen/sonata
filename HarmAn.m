function [ AnsMark ] = HarmAn( scoreMatrix, partitionNum, unmarked)
% HarmAn algorithm 找最佳路徑的方法
% input : scoreMatrix  -> 每個可能片段的分數
%         partitionNum -> 分割點數量
%         unmarked     -> 路徑可能的vertex
    if nargin < 3; unmarked = 1:partitionNum;   end

    now         = 2; 
    MARK        = 1; 
    DEL         = [];
    AnsMark(1)  = 1;
    
    while unmarked(now) < partitionNum
        if scoreMatrix(unmarked(now-1),unmarked(now))+scoreMatrix(unmarked(now),unmarked(now+1)) > scoreMatrix(unmarked(now-1),unmarked(now+1))
            MARK = [MARK, unmarked(now)];
            AnsMark(1,now) = unmarked(now);
            now = now + 1;
        else
            DEL = [DEL unmarked(now)];
            unmarked(now) = [];
        end
    end
    MARK = [MARK partitionNum];
    AnsMark(1,now) = partitionNum;


end

