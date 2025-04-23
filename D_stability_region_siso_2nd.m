function [alpha_v,beta_v,r_d,theta_s,w_H] = ...
    D_stability_region_siso_2nd(mosmax,tpmin,tsmax,csmax,trmin)
%--------------------------------------------------------------------------
% Objetivo
%   função para delimitar regiões de D-estabilidade de sistemas SISO 2a 
%   ordem a partir de valores limites de medidas transitórias especificadas 
%   nas entradas desta função.
%--------------------------------------------------------------------------
% As Entradas são:

%   "mosmax" -> valor superior limite para máximo sobressinal
%   "tpmin"  -> mínimo valor de tempo de pico
%   "tsmax"  -> máximo valor de tempo de acomodação
%   "csmax"  -> máxima porcentagem para critério de acomodação
%   "trmin"  -> mínimo valor de tempo de subida.
%--------------------------------------------------------------------------
% As Condições para entradas válidas 

%   se "mosmax" real logo 0<mosmax<=1; ou mosmax=[];
%   se "csmax" 0<csmax<1 ou csmax=[];
%   Se "trmin", "tpmin" e "tsmax" valores reais logo 0<trmin<tpmin<tsmax 
%   
%--------------------------------------------------------------------------
%   As Saídas

%	seja s a variável complexa

%   "alpha_v" -> Re{s}<=-alpha_v
%   "beta_v" -> Re{s}>=-beta_v
%   "r_d" -> abs{s}<=r_d
%   "theta_s" -> -theta_s/2<=angle{s}<=theta_s/2
%   "w_H" -> -w_H<=imag{s}<=w_H
%--------------------------------------------------------------------------
%   
    if ~isempty(tpmin)
    
        w_H = pi/tpmin;
        
    else
    
        w_H=[];
    end

    if ~isempty(mosmax)
        
        theta_s=2*atan(1/log(inv(mosmax)));
        
        if ~isempty(trmin)
            
            beta_v = (1.1*cos(theta_s/2)+...
                      0.125*cos(theta_s/2)^2+...
                      0.495*cos(theta_s)^3)/trmin;
                  
            r_d = (1.1+0.125*cos(theta_s/2)+0.495*cos(theta_s)^2)/trmin;      
        else
            
            beta_v =[];
            r_d=[];
            
        end
        if ~isempty(csmax)
       
           alpha_v = (log(csc(theta_s/2))+log(inv(csmax)))/(tsmax);
        
        else
           
           alpha_v = [];
           
       end
    else
        
        theta_s=[];
        beta_v =[];
            r_d=[];
        alpha_v = [];    
    end
        
    

end

