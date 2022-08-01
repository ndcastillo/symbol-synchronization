function [entropia letras repeticiones entropia_ll ll] = entropiaMensaje(mensaje)
    mensaje="INGENIERIA ELECTRONICA EN TELECOMUNICACIONES"
    mensaje = string(mensaje);
    letrasDivididas = split(mensaje,"");
    letrasDivididas(1) = [], letrasDivididas(end) = [];
    total = length(letrasDivididas);
    
    % Codigo ASCII
    cadena= dec2bin(char(letrasDivididas));
    
    ll = '';
    for i=1:length(cadena)    
        ll = strcat(ll,cadena(i,:));
    end
    
    total = length(ll);
    ceros = length(find(ll=='0'));
    unos = length(find(ll=='1'));
    
    entropia_ll = (ceros/total)*(log(total/ceros)/log(2)) + (unos/total)*(log(total/unos)/log(2));
    
    letras = {};
    repeticiones = [];
    probabilidades = []
    k=1;
    
    while true
        if length(letrasDivididas) == 0
            break
        else
            encontrar = find(letrasDivididas == letrasDivididas(1));
            letras{k} = letrasDivididas(1);
            repeticiones(k) = length(encontrar);
            probabilidades(k)= length(encontrar)/total;
            letrasDivididas(encontrar)=[];
        end
        k=k+1;
    end
    
    entropia = 0;
    for i=1:1:length(letras)
        entropia = entropia + repeticiones(i)*probabilidades(i)*(log(1/probabilidades(i))/log(2));
    end
    
end