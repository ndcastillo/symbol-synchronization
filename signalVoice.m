function [t F R_b B_pcm R_s B_pam quatizedSignal ...
    tb f P1 t_pam PAM tagsDec valoresCuatificacion ...
    tagsBin trama pcm_r f_pam P1_pam trama_pcm tb2 mt f_Mt Mt BPF f_sincronia ts tm numPeriodos] = signalVoice(A1,A2,A3,A4,f1,f2,f3,f4,n,L)
    A = [A1 A2 A3 A4];                    % Signal's Signal
    fm = [f1 f2 f3 f4];         % Frequency Signal
    wm = 2*pi*fm;                   % Frecquency in rad/s
    tm = 1./fm;                     % Time Period

    factor = 50;                    % Sample Factor
    frecuenciaNyquist = 2*fm;       % Nyquist Rate
    fs = factor*frecuenciaNyquist;  % Sample Frequency
    ts = 1./fs;                     % Sample Period
    rangoDinamico=5;                % Dynamic Range
    
    numMuestras = min(tm)./min(ts);
    numPeriodos = 10;
    t = 0:ts:tm*numPeriodos;

    for i=1:1:length(A) 
        x(i,:) = A(i)*cos(2*pi*fm(i).*t);
    end
    F = x(1,:);
    for i=1:1:length(A)      
        F = F + x(i,:);
    end
    
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
    
    %% PAM
    PAM = [];
    k=1;
    t_pam1 = 0:ts/numPeriodos:(d*numPeriodos-ts/numPeriodos);
    t_pam = 0:Tb_pam/(numPeriodos-1):numPeriodos*length(Fsample)*(Tb_pam/(numPeriodos-1));
    t_pam(end)=[];

    for i=1:1:length(Fsample)
        for j=1:1:numPeriodos    
            PAM(k)= Fsample(i);
            k=k+1;
        end
    end

    k=1;
    % Retention
    Fretention=reshape(Fsample,numPeriodos,[]);
    FretentionSignal = [];
    for i=1:1:length(Fretention)
        for j=1:1:numPeriodos
            FretentionSignal(k) = Fretention(1,i);
            k=k+1;
        end
    end


    %% PCM

    % Creo un Vector con los niveles de cuantificacion
    a = rangoDinamico*2/L;
    valoresCuatificacion = -5+a/2:a:5-a/2;

    % Quantizing
    quatizedSignal = FretentionSignal;
    vector = FretentionSignal;
    for i=1:1:length(FretentionSignal)
        if FretentionSignal(i) >= valoresCuatificacion(end)
            quatizedSignal(i)= valoresCuatificacion(end);
            vector(i) = L-1;
        elseif FretentionSignal(i) <= valoresCuatificacion(1)
            quatizedSignal(i)=valoresCuatificacion(1);
            vector(i) = 0;
        else
            for j=1:1:L
                if (FretentionSignal(i) > valoresCuatificacion(j) && FretentionSignal(i) < valoresCuatificacion(j) + a/2) || (FretentionSignal(i) < valoresCuatificacion(j) && FretentionSignal(i) > valoresCuatificacion(j) - a/2) 
                    quatizedSignal(i) = valoresCuatificacion(j);
                    vector(i)=j-1;
                end
            end
        end
    end

    pcm =reshape(vector,numPeriodos,[]);
    pcm_r = pcm(1,:);
    pcm_r=dec2bin(pcm_r);

    trama=[];
    trama_pcm=[];
    p=1
    k=1;

    numSamplePoints=10;
    for i=1:1:length(pcm_r)
         for j=1:1:n
             trama_pcm(p) = str2num(pcm_r(i,j));
             for d=1:1:numSamplePoints
                 trama(k) = string(pcm_r(i,j));
                 k=k+1;
             end
             p=p+1;
         end
     end

    tb=0:Tb_pcm/(numSamplePoints-1):n*numSamplePoints*length(pcm_r)*(Tb_pcm/(numSamplePoints-1));
    tb(end)=[];
    
     % Tags for Coded
            tagsDec=0:1:L-1;
            tagsBin=dec2bin(tagsDec);
            tagsBin=string(tagsBin);
            tagsBin=num2cell(tagsBin);
     % BW_pcm
            Y=fft(trama);
            P2 = abs(Y/length(trama));
            P1 = P2(1:(length(trama))/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = B_pcm*numPeriodos*(0:(length(trama)/2))/length(trama);
            
     % BW_pam
            Y_pam=fft(PAM);
            P2_pam = abs(Y_pam/length(PAM));
            P1_pam = P2_pam(1:(length(PAM))/2+1);
            P1_pam(2:end-1) = 2*P1_pam(2:end-1);
            f_pam = B_pam*numPeriodos*(0:(length(PAM)/2))/length(PAM);
            
%%----
bits_zero = zeros(1,5);
length(trama)
delay = [bits_zero trama];
trama2 = [trama bits_zero];

mt = trama2.*delay;
tb2=0:Tb_pcm/(numSamplePoints-1):n*numSamplePoints*length(pcm_r)*(Tb_pcm/(numSamplePoints-1));
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];

% mt
Y_Mt=fft(mt);
Mt2 = abs(Y_Mt/length(mt));
Mt = Mt2(1:(length(mt))/2+1);
Mt(2:end-1) = 2*Mt(2:end-1);
f_Mt = B_pcm*numPeriodos*(0:(length(mt)/2))/length(mt);

    


%% Filtro Ideal

[strong_component index_sc] = max(Mt(2:end));
% f_sincronia = round(f(index_sc+1)+(f(index_sc+2)-f(index_sc+1))/2)-20;
f_sincronia = B_pcm;
fo = find(f > f_sincronia-200,1);
vector_ff = find(f < f_sincronia+200);
ff = vector_ff(end);
isolated = f(fo:ff);
BPF = zeros(1,length(Mt));
BPF(fo:ff) = Mt(fo:ff);

end