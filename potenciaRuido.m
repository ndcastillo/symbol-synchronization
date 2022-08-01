function [x y_norm potenciaDC potenciaAC potenciaTotal entropia]=potenciaRuido(mu,desviacion)
    sigma=desviacion^2;
    x=-100:100;
    y_norm = normpdf(x,mu,sigma);
    % plot(x,y_norm)
    
    potenciaDC = mu^2;
    potenciaAC = sigma;
    potenciaTotal = potenciaDC + potenciaAC;
    entropia = log(sqrt(2*pi*exp(1))*desviacion)/log(exp(1));
end