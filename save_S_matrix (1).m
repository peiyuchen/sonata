function save_S_matrix( filename, multiProcess )
%	讀midi檔，轉成相關矩陣 左手 右手 雙手 ，休止符 ，存成.mat檔

    if nargin < 2, multiProcess = 'lowerPitch'; end  % 多音處理

     midiALL                = midi_Preprocess(filename);    % 讀ALL note資訊
    [midiRight, timeSig]    = midi_Preprocess(filename, 1); % 讀右手note資訊
     midiLeft               = midi_Preprocess(filename, 2); % 讀左手note資訊
     

    % 左右手的自相似矩陣

    if strcmp(multiProcess, 'lowerPitch') % old : 取最低音
        [ curveR, pitchInfoR ] = createCurve( midiRight, multiProcess );   % 尚未試過
        [ curveL, pitchInfoL ] = createCurve(  midiLeft, multiProcess );   % 尚未試過
        P_right = curve2pattern( pitchInfoR, curveR, timeSig );
        P_left  = curve2pattern( pitchInfoL, curveL, timeSig );
%         P_right  = pitch_Pattern_v3( midiRight, timeSig);
%         P_left   = pitch_Pattern_v3( midiLeft , timeSig);
        S.right  = similar_matrix(P_left , P_right, 'correlation', 1);
        S.left   = similar_matrix(P_left , P_left , 'correlation', 1);
        S.com    = similar_matrix(P_right, P_left , 'correlation', 1);
    elseif strcmp(multiProcess, 'allPath')    % new : 所有路徑   
        P_right  = multiPitch_Pattern_v2( midiRight, timeSig );
        P_left   = multiPitch_Pattern_v2(  midiLeft, timeSig );
        S.right  = barPath_S_matrix(P_right, P_right, 'correlation');
        S.left   = barPath_S_matrix(P_left,  P_left , 'correlation');
        S.com    = barPath_S_matrix(P_right, P_left , 'correlation');

    elseif strcmp(multiProcess, 'modifySkyline')
        [ curveT, pitchInfoT ] = createCurve( midiALL, timeSig, multiProcess, 'top'    ); % T:top
        [ curveB, pitchInfoB ] = createCurve( midiALL, timeSig, multiProcess, 'bottom' ); % B:Bottom 
        writemidi_java(pitchInfoT,['skyline_midi/' filename '_T.mid']);
        writemidi_java(pitchInfoB,['skyline_midi/' filename '_B.mid']);
        P_right  = curve2pattern( pitchInfoT, curveT, timeSig );
        P_left   = curve2pattern( pitchInfoB, curveB, timeSig );

        S.right  = similar_matrix(P_right.barPath, P_right.barPath, 'correlation', 1);
        S.left   = similar_matrix(P_left.barPath , P_left.barPath , 'correlation', 1);
        S.com    = similar_matrix(P_right.barPath, P_left.barPath , 'correlation', 1);
        midiRight= midiALL;
        midiLeft = midiALL; 
        
    end

    % 休止符
    restRight      = rest_detection( midiRight, timeSig );
    restLeft       = rest_detection(  midiLeft, timeSig );

    S.restRight    = similar_matrix( restRight, restRight, 'hamming' );
    S.restLeft     = similar_matrix(  restLeft,  restLeft, 'hamming' );
    S.restCom      = similar_matrix( restRight,  restLeft, 'hamming' );

    matrixLen = max([length(S.right), length(S.left), length(S.restRight), length(S.restLeft)]) + 1;

    S.right     ( matrixLen, matrixLen ) = 0;
    S.left      ( matrixLen, matrixLen ) = 0;
    S.com       ( matrixLen, matrixLen ) = 0;      
    S.restRight ( matrixLen, matrixLen ) = 0;
    S.restLeft  ( matrixLen, matrixLen ) = 0;
    S.restCom   ( matrixLen, matrixLen ) = 0;

    % 左右手不變自相似矩陣    
    S.right = S.right .* S.restRight;
    S.left  = S.left  .* S.restLeft ;
    S.com   = S.com   .* S.restCom  ;

    S.max2 = max( S.right, S.left ); 
    S.max3 = max( S.max2 , S.com  );

    save(['S_matrix/' multiProcess '/' filename '.mat'], 'P_right', 'P_left', 'S');

end

