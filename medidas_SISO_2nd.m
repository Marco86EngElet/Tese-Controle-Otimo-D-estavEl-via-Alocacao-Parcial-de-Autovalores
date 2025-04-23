function [mos,ts,td,tr,tp] = medidas_SISO_2nd(zeta,wn,cs)
% Funcao para obter medidas transitorias de sistemas siso de segunda ordem

%   1) Modelo SISO

%   G(s) = wn^2/(s^2+2*zeta*wn*s+wn^2);

%   2) Entradas 
%       a) zeta -> coeficiente de amortecimento
%       b) wn -> frequencia natural nao amortecida
%       c) cs -> criterio de acomodacao -> 0 < cs < 1

%   3) Saidas
%       a) mos -> maximo sobressinal
%       b) ts  -> tempo de acomodacao
%       c) td  -> tempo de atraso
%       d) tr  -> tempo de subida
%       e) tp  -> tempo de pico

%   4) Computar saidas

    mos = exp(-pi*zeta/(sqrt(1-zeta^2)));
    ts = -log(cs*sqrt(1-zeta^2))/(wn*zeta);
    td = (1.1+0.123*zeta+0.495*zeta^2)/(wn);
    tr = (1.-0.416*zeta+2.917*zeta^2)/(wn);
    tp = pi/(wn*sqrt(1-zeta^2));
end

