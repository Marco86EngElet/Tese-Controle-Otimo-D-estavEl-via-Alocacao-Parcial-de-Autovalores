%% 35) Tabela para Mostrar especificacoes transitorias e parametros das 
%       regioes transitorias

disp('Limites das Medidas Transitorias')
    Mosmax,
    
    Trmin,
    
    Trmax,
    
    Tpmin,
    
    Tsmin,
    
    Tsmax,
    
    Tdmin,
disp('------------------------------------')
disp('Parametros Calculados')
    
    rd1,
    
    rd2,

    bv1,

    bv2,

    av1,

    av2,

    av3,

disp('------------------------------------')
disp('Parametros escolhidos para D-estabilidade')
    w_H,
    
    theta_s,
    
    r_d,
    
    q_d,
    
    beta_v,
    
    alpha_v,
    
    e_P,
disp('----------------------------------')    
    
    %% 35) Tabela para Comparar Polos do sistema em malha aberta e fechada 
%   no plano das variaveis complexas.
    
    Table_Poles=table(Polos_Malha_Aberta,...
               Polos_Gcd_classico,...
               Polos_Gcd_Parcial);
    Table_Poles,
    
%% 36) Mostrar Tabelas dos Polos e Medidas Transitorias Associadas

    % 36.1)Malha Aberta
    
    Tabela_Medidas_Transitoria_Malha_Aberta,
    
    % 36.2) Projeto de alocacao otima via metodo classico
    Tabela_Medidas_Transitoria_classico,

    % 36.3) Projeto de alocacao otima via metodo parcial otimo multietapa
    Tabela_Medidas_Transitoria_Parcial,
    
%% 37) Construir tabelas para comparar normas

    %37.1) Construir coluna para normas H_2
    norma_H2=...
        [norma_H2_Godys;norma_H2_Gcdys_classico;norma_H2_Gcdys_parcial];
    
    %37.2) Construir coluna para discrepancia das normas H_2 em relacao a norma
    %   H_2 do sistema em malha aberta
    discrepancia_norma_H2_percentual=norma_H2-norma_H2_Godys;
    
    discrepancia_norma_H2_percentual=...
        100*discrepancia_norma_H2_percentual/norma_H2_Godys;
    
    %37.3) Construir coluna para normas H_infinito
    norma_Hinf=...
        [norma_Hinf_Godzs;norm_Hinf_Gcdzs_classico;norma_Hinf_Gcdzs_parcial];
    
    %37.4)Construir coluna para discrepancia das normas H_infinito em relacao 
    %   a norma H_infinito do sistema em malha aberta
    discrepancia_norma_Hinf_percentual=norma_Hinf-norma_Hinf_Godzs;
    
    discrepancia_norma_Hinf_percentual=...
        100*discrepancia_norma_Hinf_percentual/norma_Hinf_Godzs;
    
    %37.5) Construir coluna para especificar projetos de controle
    Projeto=["Malha Aberta";"Classico";"Multietapa parcial"];
    
    %37.6) Tempo de otimizacao de cada projeto
    Tempo_otimizacao=...
        [NaN('single');tempo_otimizacao_classico;tempo_otimizacao_parcial];
  
      Tabela_Normas=...
        table(Projeto,norma_H2,discrepancia_norma_H2_percentual,...
              norma_Hinf,discrepancia_norma_Hinf_percentual,...   
              Tempo_otimizacao),    