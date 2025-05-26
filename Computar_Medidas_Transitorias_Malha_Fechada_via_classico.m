%% A) Obter Coeficientes de Amortecimento, Frequencia Natural e Medidas 
%   transitorias associadas

[wn_classico,zeta_classico,p_classico]=damp(Gcdys_classico);

for i=1:size(p_classico,1)

    if zeta_classico(i)<1
        
        %A.1) Maximo Sobressinal para fracoes parciais com polos complexos
        
        Mos_Gcd_clas(i,1)=...
            exp(-pi*zeta_classico(i,1)/(sqrt(1-zeta_classico(i,1)^2)));

        %A.2) Tempo de Acomodacao para fracoes parciais com polos complexos
        
        cs = 0.05;
        Ts_Gcd_clas(i,1) =...
            -log(cs*sqrt(1-zeta_classico(i,1)^2))/(wn_classico(i,1)*zeta_classico(i,1));

        %A.3) Tempo de Atraso para fracoes parciais com polos complexos
        
        Td_Gcd_clas(i,1) =...
            (1.1+0.123*zeta_classico(i,1)+0.495*zeta_classico(i,1)^2)/(wn_classico(i,1));

        %A.4) Tempo de Subida para fracoes parciais com polos complexos
        
        Tr_Gcd_clas(i,1) =...
            (1.-0.416*zeta_classico(i,1)+2.917*zeta_classico(i,1)^2)/(wn_classico(i,1));

        %A.5) Tempo de Pico para fracoes parciais com polos complexos
        
        Tp_Gcd_clas(i,1) = pi/(wn_classico(i,1)*sqrt(1-zeta_classico(i,1)^2));
        
    else
        if zeta_classico(i)==1
            
            %A.6) Maximo Sobressinal para fracoes parciais com polo real
            
            Mos_Gcd_clas(i,1)=0;

            %A.7) Tempo de Acomodacao para fracoes parciais com polo real
            
            Ts_Gcd_clas(i,1) = 5/wn_classico(i);

            %A.8) Tempo de Atraso para fracoes parciais com polo real
            
            Td_Gcd_clas(i,1) = NaN(1);

            %A.9) Tempo de Subida para fracoes parciais com polo real
            
            Tr_Gcd_clas(i,1) = 2.2/wn_classico(i);

            %A.10) Tempo de Pico para fracoes parciais com polo real
            
            Tp_Gcd_clas(i,1) = NaN(1);
        
        else
            %A.11) Caso com polos instaveis.
            Mos_Gcd_clas(i,1)=NaN(1);
            Ts_Gcd_clas(i,1) = NaN(1);
            Td_Gcd_clas(i,1) = NaN(1);
            Tr_Gcd_clas(i,1) = NaN(1);
            Tp_Gcd_clas(i,1) = NaN(1);
        end
    end
end

%% B) Construir Tabela para Mostrar Polos e Medidas Transitorias Associados
%   do sistema LCTI-MIMO-CRPE projetado via metodo classicos com LMIs.

    %B.1) Selecionando apenas um elemento de cada par de polos complexos
    %conjugados
    
    Polos=p_classico;
    
    %B.2) Frequencia natural dos polos selecionados
    
    Wn=wn_classico;
    
    %B.3) Coeficiente de Amortecimento dos polos selecionados
    
    Zeta=zeta_classico;
    
    %B.4) Maximo Sobresinal polos selecionados
    
    Mos=Mos_Gcd_clas;
    
    if ~isempty(Mosmax)
        Mos_Ok=logical(Mos<Mosmax);    
    else
        Mos_Ok=NaN(size(Mos));
    end

    %B.5) Tempo de acomodacao polos selecionados
    
    Ts=Ts_Gcd_clas;

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
    
    %B.6) Tempo de atraso polos selecionados
    
    Td=Td_Gcd_clas;
    
    if ~isempty(Tdmin)
        Td_Ok=logical(Tdmin<=Td);
    else
        Td_Ok=NaN(size(Td));
    end
    
    %B.7) Tempo de subida polos selecionados
    
    Tr=Tr_Gcd_clas;
    
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
    
    %B.8) Tempo de pico polos selecionados
    
    Tp=Tp_Gcd_clas;
    
    if ~isempty(Tpmin)
        Tp_Ok=logical(Tpmin<=Tp);
    else
        Tp_Ok=NaN(size(Tp));
    end
    
    %B.9) Construir Tabelas relacionado Polos selecionados e medidas
    %transitorias associadas ao projeto de alocacaao otima classico
    Tabela_Medidas_Transitoria_classico=...
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
