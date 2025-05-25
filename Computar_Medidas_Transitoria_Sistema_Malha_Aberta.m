%% 9) Computar as medidas transitorias para as funcoes parciais         

    for c=1:n_x %Laco para computar medidas transitoria para cada fracao
        
        %9.1) Maximo Sobressinal para fracoes parciais com polos complexos
    
        Mos_God(c,1)=exp(-pi*zeta_God(c,1)/(sqrt(1-zeta_God(c,1)^2)));
        
        %9.2) Tempo de Acomodacao para fracoes parciais com polos complexos
        
        cs = 0.05;
        Ts_God(c,1)=...
            -log(cs*sqrt(1-zeta_God(c,1)^2))/(wn_God(c,1)*zeta_God(c,1));
    
        %9.3) Tempo de Atraso para fracoes parciais com polos complexos
        
        Td_God(c,1)=...
            (1.1+0.123*zeta_God(c,1)+0.495*zeta_God(c,1)^2)/(wn_God(c,1));
    
        %9.4) Tempo de Subida para fracoes parciais com polos complexos
        
        Tr_God(c,1)=...
            (1.-0.416*zeta_God(c,1)+2.917*zeta_God(c,1)^2)/(wn_God(c,1));
    
        %9.5) Tempo de Pico para fracoes parciais com polos complexos
        
        Tp_God(c,1)=pi/(wn_God(c,1)*sqrt(1-zeta_God(c,1)^2));
    
    end
    
%% 10) Construir Tabela para Verificar medidas transitoria associadas aos 
%     polos complexos conjugados
    
    %10.1) Polos escritos na forma de numero complexo
    
    Polos=...
        zeta_God.*wn_God+1i*wn_God.*sqrt(1-zeta_God.*zeta_God);
    
    %10.2) Frequencia natural dos polos selecionados

    Wn=wn_God;

    %10.3) Coeficiente de Amortecimento dos polos selecionados

    Zeta=zeta_God;
    
    %10.4) Maximo Sobresinal polos selecionados
    
    Mos=Mos_God;
    
    if ~isempty(Mosmax)
        Mos_Ok=logical(Mos<Mosmax);    
    else
        Mos_Ok=NaN(size(Mos));
    end
    
    %10.5) Tempo de acomodacao polos selecionados
    
    Ts=Ts_God;

    if ~isempty(Tsmin)
       if ~isempty(Tsmax)
           Ts_Ok=logical(Tsmin<=Ts<=Tsmax); 
       else
           Ts_Ok=logical(Tsmin<=Ts); 
       end
    else
       if ~isempty(Tsmax)
           Ts_Ok=logical(Ts<=Tsmax);
       else
           Ts_Ok=NaN(size(Ts));
       end
    end
    
    %10.6) Tempo de atraso polos selecionados
    
    Td=Td_God;
    
    if ~isempty(Tdmin)
        Td_Ok=logical(Tdmin<=Td);
    else
        Td_Ok=NaN(size(Td));
    end

    %10.7) Tempo de subida polos selecionados

    Tr=Tr_God;

    if ~isempty(Trmin)
        if ~isempty(Trmax)
            Tr_Ok=logical(Trmin<=Tr<=Trmax);
        else
            Tr_Ok=logical(Trmin<=Tr);
        end
    else
        if ~isemtpy(Trmax)
            Tr_Ok=logical(Tr<=Trmax);
        else
            Tr_Ok=NaN(Tr);
        end
    end
    
    %10.8) Tempo de pico polos selecionados
    
    Tp=Tp_God;

    if ~isempty(Tpmin)
        Tp_Ok=logical(Tpmin<=Tp);
    else
        Tp_Ok=NaN(size(Tp));
    end

    %10.9 Construir Tabelas relacionado Polos selecionados e medidas
    %transitorias associadas ao projeto de alocacao parcial otima
    %multietapa
    
    Tabela_Medidas_Transitoria_Malha_Aberta=...
        table(Polos,...
        Wn,...
        Zeta,...
        Mos,...
        Mos_Ok,...
        Ts,...
        Ts_Ok,...
        Td,...
        Td_Ok,...
        Tr,...
        Tr_Ok,...
        Tp,...
        Tp_Ok);