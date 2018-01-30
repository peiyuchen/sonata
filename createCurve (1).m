function [ curve, pitchInfo ] = createCurve( midiINFO, timeSig, type, isTopPitch )
%% ��midi��T�ܦ��۫ߦ��u
% input  : 1. midi��T
%          2. �縹��T
%          3. �e���ئ��u : (1) lowerPitch, (2) modifySkyline(�̰����B�̧C��)
%          4. �p�G�O modifySkyline�A�P�_�O�n�e�̰����٬O�̧C��
% 
% output : 1. �۫ߦ��u
%          2. ���u���ϥΪ�note��T

    if nargin < 3, type        = 'lowerPitch';  end
    if nargin < 4, isTopPitch  = 'top';         end
    

    innerPoint      = 100;
    
    beatNumLastSeg = (midiINFO(end,1)+midiINFO(end,2)-timeSig(end,5));
    beatNumEachBar = (timeSig(end,1)*(4/(2^timeSig(end,2))));
    totalBeat      = ceil(beatNumLastSeg/beatNumEachBar) * beatNumEachBar + timeSig(end,5);
    curveLen       = totalBeat * innerPoint;

    if strcmp(type, 'lowerPitch')

       curve = zeros(curveLen, 1);%[(eps+1) * ones(curveLen, 1)];
            
        for j = 1:size(midiINFO,1)
            Start = max(round(midiINFO(j,1)*innerPoint), 1);     % �M1��max��]�O�]��Start����O0����
            END   = round( (midiINFO(j,1) + midiINFO(j,2) ) * innerPoint );
            curve( Start : END, 1 ) = midiINFO(j,4);
        end

        pitchInfo = midiINFO;

    elseif strcmp(type, 'modifySkyline')
    
        pitchImg = zeros( max(midiINFO(:,4)), curveLen );

        for i=1:size(midiINFO,1) 
            range = max(round(midiINFO(i,1)*innerPoint),1):round((midiINFO(i,1)+midiINFO(i,2))*innerPoint);
            pitchImg(midiINFO(i,4), range ) = midiINFO(i,4);
        end

        if strcmp(isTopPitch, 'top')
            pitchImg(pitchImg == 0) =  -inf;
            curve                   =  max(pitchImg); 
            curve  (curve  ==-inf)  =  0;
        elseif strcmp(isTopPitch, 'bottom')
            pitchImg(pitchImg == 0) =  inf; 
            curve                   =  min(pitchImg); 
            curve  (curve  == inf)  =  0; 
        end
            diffPitch               =  diff(curve);
            pitchInfo               =  find(diffPitch~=0) + 1; 
            pitchInfo               =  [1, pitchInfo];
            pitchInfo(2,1:end-1)    =  diff(pitchInfo);   
            pitchInfo(2,end)        =  length(curve) + 1 - pitchInfo(1,end);
            pitchInfo(3,:)          =  1;
            pitchInfo(4,:)          =  curve(pitchInfo(1,:)); 
            pitchInfo(1:2,:)        =  pitchInfo(1:2,:) / innerPoint;
            pitchInfo(5,:)          =  120;         
            pitchInfo               =  pitchInfo(:, pitchInfo(4,:)~=0)';
            
            curve = curve';
    end

end

