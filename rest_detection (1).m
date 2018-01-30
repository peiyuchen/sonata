function [ downRest ] = rest_detection( midiInfo, timeSig )
% ヰゎ才矪瞶

%% testing code
%     clear all; close all; clc;
%     [ midiInfo, timeSig ] = midi_Preprocess('../midi/b_16_1', 2);    % 惠璶 midiInfo 1(onset beat) 4(pitch)

%% 

% initial
restLen    = ((1/2)^4); % beat
innerPoint = 100;
finalPoint =  50;
% restPos    = zeros( ceil( midiInfo(end,1) + midiInfo(end,2) ) * innerPoint , 1 );

[ pitchCurve, ~ ] = createCurve( midiInfo, timeSig );
% figure; plot( pitchCurve );

accuBar    = 0;
accuBeat   = 0;

for i = 1 : size( timeSig, 1 )    
    NoOfBeats       =   timeSig(i, 1) * ( 4 / ( 2 ^ timeSig(i, 2) ) );
    if i ~= size( timeSig, 1 )
        NoOfBar     =   ceil( ( timeSig(i+1, 5) - timeSig(i, 5) ) / NoOfBeats );
    else
        NoOfBar     =   ceil( ( midiInfo(end, 1) + midiInfo(end, 2) - timeSig(i, 5) ) / NoOfBeats );
    end

    for j = 1 : NoOfBar
        restInfo = zeros( NoOfBeats * innerPoint, 1 );
        
        barStart       = max( ( accuBeat + (j-1) * NoOfBeats ) * innerPoint, 1);
        barEnd         =      ( accuBeat +     j * NoOfBeats ) * innerPoint - 1;
        if barEnd > length(pitchCurve)
            pitchCurve = [ pitchCurve ; zeros( barEnd - length(pitchCurve), 1 ) ];
        end
        accuBar+j;
        t              = pitchCurve( barStart : barEnd )' == 0 ;
        restStart      = find( diff( [false t] ) ==  1 );
        restEnd        = find( diff( [t false] ) == -1 );
        consRest       = restEnd - restStart + 1; % 硈尿碭0       
        longRestIdx    = find( consRest > restLen * innerPoint );
        
        for idx = 1 : length( longRestIdx )
            restInfo( restStart(longRestIdx(idx) ) : restEnd(longRestIdx(idx) ), 1 ) = 1;
        end
        
        % 矪瞶竊翴计ゑ 程翴计ぶず畉翴down sample
        if NoOfBeats * innerPoint < finalPoint
            tmp = 1 / (ceil(finalPoint / (NoOfBeats * innerPoint) + 1));
            restInfo = restInfo(1 : tmp : end);
        end
        downsample                    = length(restInfo) / finalPoint;
        downRest( : , accuBar + j )   = restInfo( 1 : downsample : end, 1 );
        restPercent( accuBar + j )    = sum( downRest( : , accuBar + j ) ) / size( downRest, 1 );
        
    end
    
    accuBar  = accuBar  + NoOfBar;  
    accuBeat = accuBeat + NoOfBar * NoOfBeats;
     
end

% restSimilar = 1 - pdist2(downRest', downRest', 'hamming');


end