function [R_b Tb_pcm R_s Tb_pam] = calculateParametersVoice(ts,tm,fs,fm,t,F,numPeriodos,numMuestras)
    % Signal Voice Configuration
    ts = min(ts);
    tm = min(tm);
    d=tm;
    fs = max(fs);
    fm = max(fm);

    % Normalization
    F = F/max(F);
    F = F*rangoDinamico;

    % Square Signal
    squareSignal = zeros(1,numPeriodos);
    squareSignal(1:1)=1;
    squareSignal = repmat(squareSignal,1,numMuestras);
    F(end)=[];
    t(end)=[];

    % Sampling
    Fsample = F.*squareSignal;

    %% Computing Parameters

    % Transmission Rate 
    R_b = n*frecuenciaNyquist(1);
    disp("a) Transmission Rate: R_s= "+R_b + " bps")

    % BandWidth PCM
    Tb_pcm = 1/R_b;
    B_pcm = 1/(2*Tb_pcm);
    disp("b) BandWidth (PCM): B_pcm= "+B_pcm+" Hz")


    % 64-PAM Rate
    k = log(L)/log(2);
    R_s = R_b/k;
    disp("c) 64-PAM Rate: R_s= "+R_s + " baudios")

    % Ancho de banda PAM
    Tb_pam = 1/R_s;
    B_pam = 1/(2*Tb_pam);
    disp("d) BandWidth (PAM): B_pam = "+ B_pam + " Hz")
end