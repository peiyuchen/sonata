function [ noteBar ] = triDetection( noteBar )
%   ��tri�X��,tri�G��ӭ��ֳt�W�U����
%   ����1 : �H32�����Ŭ��ɭȡA�̫�@�ӭ����O�C
%   ����2 : ��ӭ����ʡA�������|�t�ܻ�

    channel = unique(noteBar(:,8));
    
    for ch = 1:length(channel)
        CHidx = find(noteBar(:,8)==channel(ch));
        oneTrackNote  = noteBar( CHidx, : );                               % �P�@�ӭy�D��
        dur           = round(oneTrackNote(:,2), 4);                       % ��s��32������
%         have32nd      = find(dur==0.1250);
        have32nd      = find(dur<=0.1250);
        if ~isempty(have32nd)
            t             = [diff(have32nd); 2];
            note32nd      = [];
            note32nd(:,1) = diff([0; find(t~=1)]);        
            note32nd(:,2) = have32nd([0; find(diff(have32nd)~=1)]+1);      % �s��X�� index
            note32nd(note32nd(:,1)<3,:) = [];                              % �s��3�ӥH�W          

            if ~isempty(note32nd)
                for c = 1:size(note32nd,1)
                    note32ndPitch = oneTrackNote(note32nd(c,2) : sum(note32nd(c,1:2)), 4);
                    if unique(abs(diff(note32ndPitch))) <= 3               % pitch num �t 2
                        noteBar(CHidx( note32nd(c,2)), [2 7]) = sum(oneTrackNote(note32nd(c,2) : sum(note32nd(c,1:2)), [2 7]));%noteBar(CHidx(note32nd(c,2)), [2 7]) .* (note32nd(c,1)+1); % ��tri�X��
                        noteBar(CHidx( note32nd(c,2)+1 : sum(note32nd(c,1:2)) ), :) = [];
                    end
                end
            end
        end
    end
end

