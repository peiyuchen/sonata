function [ Result ] = authenticCadence( pitch, onset, chord )
% 找 PAC
% input  : 1.和弦片段中，onset的pitch 2.和弦片段的onset 3. 屬於哪和弦
% output : 符合哪些條件

    template            = [ 0,4,7,-1;  0,4,7,10;   0,3,7,-1;    0,3,6, 9;  0,3,6,10;   0,3,6,-1 ];                      
    degree              = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 7];
    template            = ceil(chord/12);
    root                = ~(ceil(mod(chord,12)/12))*12 + mod(chord,12);
    root(template==0)   = 0;
    
    % 把root變成一維，紀錄每個root為哪個bar (chordBar)
    pitchOneD  = reshape(pitch'   ,1,[]);
    onsetOneD  = reshape(onset'   ,1,[]);
    tempOneD   = reshape(template',1,[]);
    rootOneD   = reshape(root'    ,1,[]);
    rootBar    = reshape(repmat((1:size(root))',1,size(root,2))',1,[]);
      
    Null       = rootOneD==0;
    rootBar   (Null) = [];
    onsetOneD (Null) = [];
    rootOneD  (Null) = [];
    tempOneD  (Null) = [];
    pitchOneD(cellfun('size', pitchOneD, 1)==0) = [];


    %% init
    LEN         = length(rootOneD);
    V2I         = zeros(LEN,1);
    bigV        = zeros(LEN,1);
    VrootPos    = zeros(LEN,1);
    IrootPos    = zeros(LEN,1); 
    IhighPos    = zeros(LEN,1);
    Idownbeat   = zeros(LEN,1);
    
    %% V -> I  
    rootDegree      = degree(rootOneD);
    pitchClassDiff  = diff(rootDegree);
    pitchClassDiff(pitchClassDiff>0) = pitchClassDiff(pitchClassDiff>0) - 7;
    V2I_idx         = find(pitchClassDiff==-4);      
    V2I(V2I_idx)    = 1;   
    
    %% V : 大和弦
    bigV_idx        = V2I_idx(tempOneD(V2I_idx)<3);  
    bigV(bigV_idx)  = 1;	
   
    %% 以下為3個條件 PAC,不符合任一為IAC
    idx = bigV_idx
    
    for i=1:length(idx)
        VlowPitch (i)       = mod(pitchOneD{idx(i)  }(1)  ,12) + 1; 
        IlowPitch (i)       = mod(pitchOneD{idx(i)+1}(1)  ,12) + 1;
        IhighPitch(i)       = mod(pitchOneD{idx(i)+1}(end),12) + 1;
        
        VrootPos (idx(i))   = (VlowPitch(i)  == rootOneD(idx(i)  )); % (1) V 原位
        IrootPos (idx(i))   = (IlowPitch(i)  == rootOneD(idx(i)+1)); % (1) I 原位
        IhighPos (idx(i))   = (IhighPitch(i) == rootOneD(idx(i)+1)); % (2) I 高八度   
        Idownbeat(idx(i))   = (onsetOneD(idx(i)+1) == 0);            % (3) I 強拍 

        TRUE_IrootPos(i)    = (IrootPos(idx(i)) && IrootPos(idx(i)) && VrootPos (idx(i))) * idx(i);
        TRUE_ALL(i)         = (VrootPos(idx(i)) && IrootPos(idx(i)) && IhighPos(idx(i)) && Idownbeat(idx(i))) * idx(i);
    end
    
    TRUE_IrootPos(TRUE_IrootPos==0) = [];
    TRUE_ALL(TRUE_ALL==0)           = [];
    
    % save to result
    Result.rootBar    = rootBar';
    Result.tempOneD   = tempOneD';
    Result.rootOneD   = rootOneD';
    Result.rootDegree = rootDegree';
    Result.V2I        = V2I;
    Result.bigV       = bigV;
    Result.VrootPos   = VrootPos;
    Result.IrootPos   = IrootPos;
    Result.IhighPos   = IhighPos;
    Result.Idownbeat  = Idownbeat;
    
    Result            = struct2table(Result);
    
    
    PAC_bar = [Result.rootBar(TRUE_ALL,1)      Result.rootBar(TRUE_ALL+1,1)     ]
    IAC_bar = [Result.rootBar(TRUE_IrootPos,1) Result.rootBar(TRUE_IrootPos+1,1)]
    
    
end

