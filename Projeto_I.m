clc, close, clear

format SHORT
%% 1) Definir quantidade de entradas e saidas.

    % 1.1) Quantidade de estados
    n_x=10;
    
    % 1.2) Quantidade de saidas nao sem influencia direta de disturbios
    n_y=1;
    
    % 1.3) Quantidade de saidas com influencia direta de disturbios
    n_z=1;

    % 1.4) Quantidade de de disturbios
    n_d=1;    

    % 1.5) Quantidade de acoes de controle
    n_u=2;

    % 1.6) Quantidade autovalores a serem modificados.
    n_ola=n_u;
 
%% 2) Inserir a quantidade de polos reais $n_r$ e par de polos complexos 
%   conjugados $n_{pcc}$.

    %2.1) Quantidade de polos reais
    n_r=0;
    
    %2.2) Quantidade de pares de polos complexos
    n_pcc=5;

%% 3) Delimitar limites para medidas transitoria associadas:

    %3.1) Sobressinal maximo.
    
    Mosmax=0.2;
    
    %3.2) Tempo de subida minimo e maximo.
    
    Trmin=1.8;
    Trmax=6;
    
    %3.3) Tempo de pico minimo.
    
    Tpmin=2.5;
    
    %3.4) Tempo de acomodacao maximo e minimo.
    
    csmax=0.02;
    Tsmax=8.8;
    Tsmin=5;
    
    %12.5) Tempo de atraso minimo.
    
    Tdmin=[];

%% 4) Computar os parametros das LMIs de D-estabilidade:

    Computar_parametros_D_estabilidade
    
%% 5) Computar as funcoes de transferencia

    %5.1) Entradas para coeficiente de amortecimentos das fracoes parciais
    
    zeta_c(1,1) = 0.1;  % Funcao parcial 1
    zeta_c(2,1) = 0.12; % Funcao parcial 2
    zeta_c(3,1) = 0.15; % Funcao parcial 3
    zeta_c(4,1) = 0.2;  % Funcao parcial 4
    zeta_c(5,1) = 0.3;  % Funcao parcial 5

    %5.2) Entradas para frequencia natural das fracoes parciais
    
    w_c(1,1) = 0.1;  % Funcao parcial 1
    w_c(2,1) = 0.12; % Funcao parcial 2
    w_c(3,1) = 0.15; % Funcao parcial 3
    w_c(4,1) = 0.2;  % Funcao parcial 4
    w_c(5,1) = 0.3;  % Funcao parcial 5

    %5.3) Entradas para Ganho das fracoes parciais
    
    KYc=ones(5,1); %Ganho das fracoes parciais associadas a "Gody(s)" 
    KZc=ones(5,1); %Ganho das fracoes parciais associadas a "Godz(s)"
    KZc(3,1)=-1; %Correcao de Ganho de uma das fracoes parciais associadas a "Godz(s)"

    %5.4) Escrita das funcoes de transferencia "Gody(s)" e "Godz(s)" 

        %5.4.1) Criar entrada para funcoes
    
        Godys=tf(0,1); % Gody(s)=0;
        Godzs=tf(0.01,1); % Godz(s)=0.01;
        
        %5.4.2) Estrutura de laco para somar fracoes parciais 
        
        for c=1:5 % Contador de laco
            num = w_c(c,1)^2; % Numerator da fracao parcial
            den = [1, 2*zeta_c(c,1)*w_c(c,1), w_c(c,1)^2]; % Denominator da fracao parcial 
            G_c = tf(num, den); % Fracao parcial sem ganho
            Godys=Godys+KYc(c,1)*G_c; % Adicionar fracao parcial em "Gody(s)"
            Godzs=Godzs+KZc(c,1)*G_c; % Adicionar fracao parcial em "Godz(s)"
        end

%% 6) Representacao canonica das funcoes de transferencia do sistema em malha
%    Aberta
    
    %6.1) Computa representcao canonica controlavel de Gody(s) 
    
    Canon_Godys = ss(Godys);
    
    %6.2) Computa representacao canonica controlavel de Godz(s)
    
    Canon_Godzs = ss(Godzs);
    
    %6.3) Computar as matrizes do sistema LCTI-MIMO-CRPE em malha aberta
    
    A_x=Canon_Godys.A;
    B_d=Canon_Godys.B;
    C_y=Canon_Godys.C;
    E_y=Canon_Godys.D;
    C_z=Canon_Godzs.C;
    E_z=Canon_Godzs.D; 
    
%% 7) Computar Normas e Polos do sistema em Malha Aberta
    
    Computar_Normas_e_Polos_Malha_Aberta

%% 8) Computar Frequencia Natural e Coeficiente de amortecimento

    [wn_God,zeta_God,p_God]=damp(Canon_Godys);
    
%% 9) Computar Medidas transitorias do sistema em Malha Aberta

    Computar_Medidas_Transitoria_Sistema_Malha_Aberta
    
%% 10) Computar respostas de Gody(s) ao disturbio impulso unitario

    tempo_simulacao_impulso=0:0.01:14;
    
    Computar_Dados_Simulacao_impulso
    
%% 11) Especificar as Matrizes associadas ao controle 'B_u', 'D_y', 'D_z'.
    
    %11.1) B_u
    
    B_u=zeros(n_x,n_u);
    B_u(2,1)=1;
    B_u(3,1)=1;
    B_u(2,2)=-1;
    B_u(3,2)=-1;
    
    %11.2) D_y e D_z
    
    D_y=0.001*ones(n_y,n_u);
    D_z=D_y;
    
%% 12) Configuracoes no yalmip

    %12.1) Configuracoes para otimizacao semidefinida no yalmip
    
    Configuracoes_SDP =...
         sdpsettings('verbose',1,'solver','lmilab','debug',1);

    %12.2) Pesos da funcao custo de otimizacao 
    
    c_H2=1;
    c_Hinf=2;

%% 13) Computar quantidade de estados nao alocaveis e checar controlabilidade total
    
    Computar_quantidade_de_estados_nao_alocaveis

%% 14)  Criar matrizes e variaveis de decisao para LMIs 

X  = sdpvar(n_x,n_x,'symmetric');
W  = sdpvar(n_u,n_x,'full');
gamma = sdpvar(1,1,'symmetric'); 
Z = sdpvar(n_y,n_y,'symmetric');
rho = sdpvar(1,1,'symmetric');   
        
%% 15) Computar as LMIs para controle otimo misto D-estavel

LMIs_controle_misto_D_estavel

%% 16) Executar Algoritmo 3

tic
optimize(set_LMIs_classico,...
            c_H2*rho+c_Hinf*gamma,...
            Configuracoes_SDP);
tempo_otimizacao_classico=toc;

%% 17) Extrair variaveis otimas do Algoritmo 3

optimal_W = value(W);
optimal_X = value(X);

%% 18) Computar Matriz de retroalimentacao via Algoritmo 3

Kpf_classico=optimal_W/optimal_X;

%% 19) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 3

    % 19.1) Espaco de Estados para Gcdy(s) 
    
    Gcdys_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_y+D_y*Kpf_classico,E_y);

    % 19.2) Espaco de Estados para Gcdz(s)  
    
    Gcdzs_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_z+D_z*Kpf_classico,E_z);

%% 20) Computar medidas transitorias e tabela para sistema em malha fechada projetado via algoritmo 3 

    Computar_Medidas_Transitorias_Malha_Fechada_via_classico
    
%% 21) Computar Normas H_2 e H_infinito para sistema projetado via 
%   Algoritmo 3

    %21.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_classico=norm(Gcdys_classico,2);
    
    %21.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norm_Hinf_Gcdzs_classico=norm(Gcdzs_classico,'inf');
    
%% 22) Computar Os polos do sistema projetado via Algoritmo 3

    Polos_Gcd_classico=cplxpair(pole(Gcdys_classico));
    
%% 23) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_classico
    
%% 24) Executar Algoritmo 7 multietapa

    %24.1) Iniciar matriz de retroalimentacao
    
    Kpf_parcial=zeros(n_u,n_x);
    
    %24.2) Entrada Q_p  
    
    Q_p = [1 1;1i -1i];

    %24.3) Iniciar contagem de tempo para otimizacao multietapa via 
    %       Algoritmo 7
    tic
    
    %24.4) Laco para realizar otimizacao multietapa
    for j = 1:n_pcc % contagem de par de polos selecionado
        
        %24.4.1) Obter Matrizes dos Autovalores e Autovetores
        
        [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
            eig(A_x + B_u * Kpf_parcial);

        %24.4.2) Escolher o j-enesimo par de polos contidos em "Open_Loop_Poles"   
        
        pole1 = Polos_Malha_Aberta(2 * j - 1, 1);
        pole2 = Polos_Malha_Aberta(2 * j, 1);

        %24.4.3) Computar as diferencas entre os autovalores e os polos
        
        dif_pole1 = abs(Eigenvalues - pole1);
        dif_pole2 = abs(Eigenvalues - pole2);

        %24.4.4) Encontrar o indice dos elementos que tem a menor diferenca
        % entre os autovalores e os polos selecionados para alocacao
        
        [~, minIndex_1] = min(dif_pole1(:));
        [~, minIndex_2] = min(dif_pole2(:));

        %24.4.5) Converte o indices lineares "minIndex_1" e "minIndex_2"  
        % do menor valor para subscritos (linha e coluna, contidos em 
        % "row_1", "row_2", "col_1", "col_2") da matriz "Eigenvalues".
        
        [row_1, col_1] = ind2sub(size(Eigenvalues), minIndex_1);
        [row_2, col_2] = ind2sub(size(Eigenvalues), minIndex_2);

        %24.4.6) Ordenar os intervalos de linhas e colunas contidos em 
        % "row_poles" e "col_poles".   
        
        row_poles = min(row_1, row_2):max(row_1, row_2);
        col_poles = min(col_1, col_2):max(col_1, col_2);

        %24.4.7) Encontrar a submatrices "Lambda_p" e "L_p"
        
        Lambda_p = Eigenvalues(row_poles, col_poles);
        L_p = Left_Eigenvectors(1:end, col_poles);

        %24.5) Computar controle otimo misto D-estavel com alocacao parcial
        Controle_otimo_misto_D_estavel_alocacao_parcial
    
            %24.4.10.11) Executar Algoritmo 7 Multietapa   
            if ~isempty(c_H2) % Se houver otimizacao para norma H2
                 if ~isempty(c_Hinf) % Se houver otimizaco para norma Hinf
                     
                    optimize(set_LMIs_parcial,...
                        c_H2*rho+c_Hinf*gamma,Configuracoes_SDP); 
                    
                 else % Se nao houver otimizacao para norma Hinf
                     
                     optimize(set_LMIs_parcial,...
                         c_H2*rho,Configuracoes_SDP);
                     
                 end
            else % Se nao houver otimizacao para norma H2
                if ~isempty(c_Hinf) % Se houver otimizacao para norma Hinf
                    
                    optimize(set_LMIs_parcial,...
                        c_Hinf*gamma,Configuracoes_SDP); 
                    
                else % Se nao houver otimizacao para norma Hinf
                    
                     optimize(set_LMIs_parcial,[],Configuracoes_SDP);
                     
                end 
            end

            %24.4.10.12) Extrair variaveis de decisao para computar K_pf 
            
            optimal_Tilde_W = value(Tilde_W);
            optimal_Tilde_X = value(Tilde_X);

            %24.4.10.13) Computar e Atualizar K_pf
            
            KD_parcial=optimal_Tilde_W*inv(optimal_Tilde_X);
            Kpf_parcial = Kpf_parcial + KD_parcial * (Q_p * L_p') / 2;


    end %finalizar Algoritmo 7 multietapa
    
    % 24.5) Finalizar tempo de otimizacao multietapa
    
    tempo_otimizacao_parcial = toc;
    
%% 25) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 7 multietapa

    % 25.1) Espaco de estados para Gcdy(s) via algoritmo 7 multietapa
    
    Gcdys_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);

    % 25.2) Espaco de estados para Gcdz(s) via algoritmo 7 multietapa
    
    Gcdzs_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z); 

%% 26) Computar Normas H_2 e H_infinito para sistema projetado via Algoritmo 7

    %26.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_parcial=norm(Gcdys_parcial,2);
    
    %26.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norma_Hinf_Gcdzs_parcial=norm(Gcdzs_parcial,'inf');
    
%% 27) Computar Os polos do sistema via teorema 7

    Polos_Gcd_Parcial=cplxpair(pole(Gcdys_parcial));
    
%% 28) Obter Coeficientes de Amortecimento, Frequencia Natural e Medidas 
%   transitorias associadas

[wn_theo7,zeta_theo7,p_theo7]=damp(Gcdys_parcial);
cont=1;

%% 29) Computar medidas transitorias do sistema com retroalimentação projetada via algoritmo 5

Computar_Medidas_Transitorias_Malha_Fechada_via_parcial

%% 30) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_parcial
    
%% 31) Computar Tabelas
x
    Computar_Tabelas
          
%% 32) Construir graficos de Bode para sistema LCTI-MIMO em malha aberta e 
%   fechada.

    figure 
    sigmaplot(Gcdzs_classico,Gcdzs_parcial,Canon_Godzs)
    title('Valor Singular G_d_z(s)')
    legend('Gcdz(s) Classico','Gcdz(s) Parcial',...
           'Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines
    
%% 33) Construir graficos de resposta transitoria de g(t)^2 ao disturbio  
%   impulso unitario do sistema LCTI-MIMO em malha aberta e fechada.

    figure
    subplot(121)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico,...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial,...
         t_Canon_Godys_impulse,y2_Canon_Godys_impulse);
    xlim([0 14])
    title('Malha Aberta e Fechada')
    xlabel('t')
    ylabel('y^2(t)')
    legend('Gcdy(s) classico','Gcdy(s) parcial','Gody(s)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico,...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial,...
         t_Canon_Godys_impulse,y2_Canon_Godys_impulse);
    xlim([0 5])
    ylim([0 0.008])
    title('Classico Vs Parcial')
    xlabel('t')
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Resposta de y(t) ao impulso')
    
%% 34) Construir grafico de u(t) sistema LCTI-MIMO em malha aberta e 
%   fechada em resposta a disturbio impulso unitario.     
    
    figure
    subplot(121)
    semilogy(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,1),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,1));
    xlim([0,11]) 
    title('u_1(t)')
    xlabel('t')
    ylabel('u^2_1(t)')
    legend('G_c_d_y(s) classico','G_c_d_y(s) parcial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    semilogy(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,2),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,2));
    xlim([0,11]) 
    title('u_2(t)')
    xlabel('t')
    ylabel('u^2_2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Resposta de u(t) ao impulso')
    
%% 35) Grafico para Localizacao dos Polos

    %35.1) Autovalores (Polos) a serem inseridos no grafico
    autovalores=[ Polos_Malha_Aberta',...
                  Polos_Gcd_classico',...  
                  Polos_Gcd_Parcial'];

    %35.2) Computando curvas de fronteira das regioes    
    [  x_horizontal_strip,y_inferior_horizontal_strip,...
        y_superior_horizontal_strip,...
        x_sector,y_inferior_sector,y_superior_sector,...
        x_parabola,y_inferior_parabola,y_superior_parabola,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = linhas_D_regioes(...
        autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

    %35.3) Grafico            
    figure
    subplot(121)
    plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
         real(Polos_Gcd_classico),...
         imag(Polos_Gcd_classico),'+r',...
         real(Polos_Gcd_Parcial),...
         imag(Polos_Gcd_Parcial),'xb',...
         x_alphav,y_alphav,'-k',...
         x_betav,y_betav,'-k',...
         x_sector,y_inferior_sector,'-k',...
         x_sector,y_superior_sector,'-k',...
         x_disk,y_disk,'-k',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-k',...
         x_horizontal_strip,y_superior_horizontal_strip,'-k',...
         'LineWidth',2,'MarkerSize',8);
   legend('Aberta','Livre','Parcial','Fronteira','FontSize',10) 
   grid on
   
   subplot(122)
   plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
         real(Polos_Gcd_classico),...
         imag(Polos_Gcd_classico),'+r',...
         real(Polos_Gcd_Parcial),...
         imag(Polos_Gcd_Parcial),'xb',...
         x_alphav,y_alphav,'-k',...
         x_betav,y_betav,'-k',...
         x_sector,y_inferior_sector,'-k',...
         x_sector,y_superior_sector,'-k',...
         x_disk,y_disk,'-k',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-k',...
         x_horizontal_strip,y_superior_horizontal_strip,'-k',...
         'LineWidth',2,'MarkerSize',8);
   ylim([-0.6 0.6])
   xlim([-0.71 0])
   grid on
   sgtitle('Polos e localizacao Plano Complexo','FontSize',12)
