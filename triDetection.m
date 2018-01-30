function [ noteBar ] = triDetection( noteBar )
%   把tri合併,tri：兩個音快速上下跳動
%   條件1 : 以32分音符為時值，最後一個音不是。
%   條件2 : 兩個音跳動，音高不會差很遠

    channel = unique(noteBar(:,8));
    
    for ch = 1:length(channel)
        CHidx = find(noteBar(:,8)==channel(ch));
        oneTrackNote  = noteBar( CHidx, : );                               % 同一個軌道找
        dur           = round(oneTrackNote(:,2), 4);                       % 找連續的32分音符
%         have32nd      = find(dur==0.1250);
        have32nd      = find(dur<=0.1250);
        if ~isempty(have32nd)
            t             = [diff(have32nd); 2];
            note32nd      = [];
            note32nd(:,1) = diff([0; find(t~=1)]);        
            note32nd(:,2) = have32nd([0; find(diff(have32nd)~=1)]+1);      % 連續幾個 index
            note32nd(note32nd(:,1)<3,:) = [];                              % 連續3個以上          

            if ~isempty(note32nd)
                for c = 1:size(note32nd,1)
                    note32ndPitch = oneTrackNote(note32nd(c,2) : sum(note32nd(c,1:2)), 4);
                    if unique(abs(diff(note32ndPitch))) <= 3               % pitch num 差 2
                        noteBar(CHidx( note32nd(c,2)), [2 7]) = sum(oneTrackNote(note32nd(c,2) : sum(note32nd(c,1:2)), [2 7]));%noteBar(CHidx(note32nd(c,2)), [2 7]) .* (note32nd(c,1)+1); % 把tri合併
                        noteBar(CHidx( note32nd(c,2)+1 : sum(note32nd(c,1:2)) ), :) = [];
                    end
                end
            end
        end
    end
end

