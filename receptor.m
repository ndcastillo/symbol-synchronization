function [t_clock clock samples t m_estimada] = receptor(f_sincronia,lenTrama,n,trama,valoresCuatificacion,ts,tm,numPeriodos)
%% Demodulacion
B = f_sincronia;
Tb = 1/(2*B);
numSamplePoints=10;
yt = [];
a= 'Hola';
f_sincronia
lenTrama

%% data clock
vector = [ones(1,numSamplePoints/2) zeros(1,(numSamplePoints/2))]';
vector = [vector vector];
data = vector;
for i=1:1:(lenTrama/(10*2)-1)
    data = [data vector];
end

clock = reshape(data,1,[]);
clock(1)=[]
t_clock = 0:Tb/9:n*numSamplePoints*(100)*(Tb/(numSamplePoints-1));
t_clock(end)=[];t_clock(end)=[];

%% sampling
trama(end)=[]

samples = trama.*clock;

k = 1;
for i=1:1:length(samples)
   if mod(i,10) == 1
       yt(k) = samples(i)
       k = k + 1
   end
end

 % Detector de Simbolo

simbolos = string(reshape(yt,n,[])');
dim_simbolos = size(simbolos);
vector_sim =[]
volts_simbol = []

for i=1:dim_simbolos(1)
    sis=[]
    for j=1:dim_simbolos(2)
        sis = strcat(sis,simbolos(i,j));
    end
    vector_sim(i) = sis;
    volts_simbol(i)=valoresCuatificacion(bin2dec(sis)+1);
end

m_estimada=[]
k=1
for i=1:1:length(volts_simbol)
     for d=1:1:numSamplePoints
         m_estimada(k) = volts_simbol(i);
         k=k+1;
     end
end
 
t = 0:ts:tm*numPeriodos;
t(end)=[]
% 
% figure
% plot(t,m_estimada, 'LineWidth',1.5);
%     yticks(valoresCuatificacion)
%     yticklabels(tagsBin)
%     style = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',style,'fontsize',8)
%     set(gca,'XTickLabelMode','auto')
%     title('Coded Signal')
%     ylabel('Levels of Voltage [V]')
%     xlabel('t[s]')
%     grid on;
end