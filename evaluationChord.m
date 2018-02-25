clear; clc; close all;

GTFile  = 'trans_chord_k545';
evaFile = 'testBoth';
timeSig = 4;    unit = 1;
[number, text,  GTdata] = xlsread(['chordGT/' GTFile '.xlsx']);
[number, text, EVAdata] = xlsread(['chordEva/' evaFile '.xlsx']);


[GTarray ,GTname ] = toArray( GTdata, unit, timeSig);
[EVAarray,EVAname] = toArray(EVAdata, unit, timeSig);
Ans = (GTarray-EVAarray==0);
ans1 = sum(sum(GTarray-EVAarray==0))/(size(GTarray,1)*size(GTarray,2));

[GTarray ] = toArray( GTdata, unit/2, timeSig);
[EVAarray] = toArray(EVAdata, unit/2, timeSig);

ans2 = sum(sum(GTarray-EVAarray==0))/(size(GTarray,1)*size(GTarray,2));

function [chordNo, chordName] = toArray(data, unit, timeSig)
    barNum      = data{end,1};
    chordNo     = zeros(barNum, timeSig/unit);
    chordName   = repmat({'-'}, barNum, timeSig/unit);%cell (barNum, timeSig/unit);
    
    for i=2:length(data)
        on = floor(data{i,2}/unit)+1; 
        if i~=length(data)
            off = ceil(data{i+1,2}/unit);
            if data{i+1,2}==0; off = 4/unit; end
        else
            off = 4/unit;
        end
        chordNo  (data{i,1}, on:off) = data{i,5};
        chordName{data{i,1}, on}     = data{i,4};
    end
end