function [n L] = calculateSNR(snrQuatizationdB)
    %***** SNR (QUATIZED - dB) ***** 
    snrQuatization = 10^(snrQuatizationdB/10);
    L = sqrt(snrQuatization/3);     % Levels

    if L==0 || L < 0
        disp('Fuera del Rango Establecido')
    elseif L > 0 && L<=1
        L=1;
    elseif L>1 && L<=2
        L=2;
    elseif L>2 && L<=4
        L=4;
    elseif L>4 && L<=8
        L=8;
    elseif L>8 && L<=16
        L=16;
    elseif L>16 && L<=32
        L=32;
    elseif L>32 && L<=64
        L=64;
    elseif L>64 && L<=128
        L=128;
    elseif L>128 && L<=256
        L=256;
    end

    n = log(L)/log(2);
end

