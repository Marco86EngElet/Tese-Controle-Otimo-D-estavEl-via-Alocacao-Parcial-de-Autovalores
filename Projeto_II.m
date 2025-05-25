clc, close all, clear all

format SHORT

%% 1) Definir quantidade de entradas e saidas.

    % 1.1) Quantidade de estados
    n_x=12;
    
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
    n_pcc=6;

%% 3) Delimitar limites para medidas transitoria associadas:

    %3.1) Sobressinal maximo.
    Mosmax=0.06;
    
    %3.2) Tempo de subida minimo e maximo.
    
    Trmin=0.5;
    Trmax=[];
    
    %3.3) Tempo de pico minimo.
    
    Tpmin=[];
    
    %3.4) Tempo de acomodacao maximo e minimo.
    
    csmax=0.02;
    Tsmax=4.71;
    Tsmin=[];
    
    %12.5) Tempo de atraso minimo.
    
    Tdmin=[];

%% 4) Computar os parametros das LMIs de D-estabilidade:

    Computar_parametros_D_estabilidade
    
%% 5) Computar as funcoes de transferencia

    %5.1) Entradas para coeficiente de amortecimentos das fracoes parciais
    
    zeta_c(1,1) = 0.9^(1/30);  % Funcao parcial 1
    zeta_c(2,1) = 0.9^(2/30); % Funcao parcial 2
    zeta_c(3,1) = 0.9^(3/30); % Funcao parcial 3
    zeta_c(4,1) = 0.9^(4/30);  % Funcao parcial 4
    zeta_c(5,1) = 0.9^(5/30);  % Funcao parcial 5
    
    zeta_c(6,1)=0.2; %Funcao parcial a ser alocada
    
    %5.2) Entradas para frequencia natural das fracoes parciais
    
    wn_c(1,1) = exp(-2+1/30);  % Funcao parcial 1
    wn_c(2,1) = exp(-2+2/30); % Funcao parcial 2
    wn_c(3,1) = exp(-2+3/30); % Funcao parcial 3
    wn_c(4,1) = exp(-2+4/30);  % Funcao parcial 4
    wn_c(5,1) = exp(-2+5/30);  % Funcao parcial 5

    wn_c(6,1) = 0.1;  % Funcao parcial 5

    %5.3) Entradas para Ganho das fracoes parciais
    
    KYc=0.5*ones(6,1); %Ganho das fracoes parciais associadas a "Gody(s)" 
    KZc=0.5*ones(6,1); %Ganho das fracoes parciais associadas a "Godz(s)"
    
    KYc(6,1)=10; %Correcao de Ganho de uma das fracoes parciais associadas a "Godz(s)"
    KZc(6,1)=10; %Correcao de Ganho de uma das fracoes parciais associadas a "Godz(s)"

    %5.4) Escrita das funcoes de transferencia "Gody(s)" e "Godz(s)" 

        %5.4.1) Criar entrada para funcoes
    
        Godys=tf(0,1); % Gody(s)=0;
        Godzs=tf(1,1); % Godz(s)=0.01;
        GNAs=tf(0,1);%GNA(s)=0;
        %5.4.2) Estrutura de laco para somar fracoes parciais 
        
        for c=1:6 % Contador de laco
            num = wn_c(c,1)^2; % Numerator da fracao parcial
            den = [1, 2*zeta_c(c,1)*wn_c(c,1), wn_c(c,1)^2]; % Denominator da fracao parcial 
            G_c = tf(num, den); % Fracao parcial sem ganho
            Godys=Godys+KYc(c,1)*G_c; % Adicionar fracao parcial em "Gody(s)"
            Godzs=Godzs+KZc(c,1)*G_c; % Adicionar fracao parcial em "Godz(s)"
            if c<=5
                GNAs=GNAs+KYc(c,1)*G_c;
            else
                GAs=KYc(c,1)*G_c;
            end
        end

%% 6) Representacao canonica das funcoes de transferencia do sistema em malha
%    Aberta
    
    %6.1) Computa representacao canonica controlavel de Gody(s) 
    
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
    
    %13.1) B_u
    
    B_u=zeros(n_x,n_u);
    B_u(2,1)=1;
    B_u(3,1)=1;
    B_u(2,2)=-1;
    B_u(3,2)=-1;
    
    %13.2) D_y e D_z
    
    D_y=0.001*ones(n_y,n_u);
    D_z=D_y;
    
%% 14) Configuracoes no yalmip

    %14.1) Configuracoes para otimizacao semidefinida no yalmip
    
    Configuracoes_SDP =...
         sdpsettings('verbose',1,'solver','lmilab','debug',1);

    %14.2) Pesos da funcao custo de otimizacao 
    
    c_H2=2;
    c_Hinf=1;

%% 15) Computar matriz de retroalimentacao com alocacao livre

    Computar_quantidade_de_estados_nao_alocaveis

%% 16)  Criar matrizes e variaveis de decisao para LMIs 

X  = sdpvar(n_x,n_x,'symmetric');
W  = sdpvar(n_u,n_x,'full');
gamma = sdpvar(1,1,'symmetric'); 
Z = sdpvar(n_y,n_y,'symmetric');
rho = sdpvar(1,1,'symmetric');   
        
%% 17) Computar as LMIs para controle otimo misto D-estavel

LMIs_controle_misto_D_estavel

%% 20) Executar Algoritmo 3

tic
optimize(set_LMIs_classico,...
            c_H2*rho+c_Hinf*gamma,...
            Configuracoes_SDP);
tempo_otimizacao_classico=toc;

%% 21) Extrair variaveis otimas do Algoritmo 3

optimal_W = value(W);
optimal_X = value(X);

%% 22) Computar Matriz de retroalimentacao via Algoritmo 3

Kpf_classico=optimal_W/optimal_X;

%% 23) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 3

    % 23.1) Espaco de Estados para Gcdy(s) 
    
    Gcdys_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_y+D_y*Kpf_classico,E_y);

    % 23.2) Espaco de Estados para Gcdz(s)  
    
    Gcdzs_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_z+D_z*Kpf_classico,E_z);

%% 23) Computar medidas transitorias e tabela para sistema em malha fechada projetado via algoritmo 3 

    Computar_Medidas_Transitorias_Malha_Fechada_via_classico
    
%% 26) Computar Normas H_2 e H_infinito para sistema projetado via 
%   Algoritmo 3

    %26.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_classico=norm(Gcdys_classico,2);
    
    %26.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norm_Hinf_Gcdzs_classico=norm(Gcdzs_classico,'inf');
    
%% 27) Computar Os polos do sistema projetado via Algoritmo 3

    Polos_Gcd_classico=cplxpair(pole(Gcdys_classico));

%% 27) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_livre
    
%% 29) Executar Algoritmo 7 multietapa

    %29.1) Iniciar matriz de retroalimentacao 
    
    Kpf_parcial=zeros(n_u,n_x);

    %29.2) Entrada Q_p
    
    Q_p = [1 1;1i -1i];

    %29.3) Obter Matrizes dos Autovalores e Autovetores
        
    [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
        eig(A_x + B_u * Kpf_parcial);
        
    %29.4) Encontrar as matrizes para alocacao
    [row_poles,col_poles]=find(...
        real(Eigenvalues)<-0.019 & real(Eigenvalues)>-0.021)
    
    Lambda_p = Eigenvalues(row_poles, col_poles);
    L_p = Left_Eigenvectors(1:end, col_poles);    
        
    %29.5) Computar controle otimo misto D-estavel com alocacao parcial 
    Controle_otimo_misto_D_estavel_alocacao_parcial
    
        %29.18) Iniciar contagem de tempo para otimizacao multietapa via 
        %       Algoritmo 7
        tic
    
        %29.19) Executar Algoritmo 7 Multietapa   
        if ~isempty(c_H2) % Se houver otimizacao para norma H2
             if ~isempty(c_Hinf) % Se houver otimizacao para norma Hinf

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

        %29.20) Extrair variaveis de decisao para computar K_pf 

        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);

        %29.21) Computar e Atualizar K_pf

        KD_parcial=optimal_Tilde_W*inv(optimal_Tilde_X);
        Kpf_parcial = Kpf_parcial + KD_parcial * (Q_p * L_p') / 2;

        
    % 29.22) Finalizar tempo de otimizacao multietapa
    
    tempo_otimizacao_parcial = toc;
    
%% 30) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 7 multietapa

    % 30.1) Espaco de estados para Gcdy(s) via algoritmo 7 multietapa
    
    Gcdys_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);

    % 30.2) Espaco de estados para Gcdz(s) via algoritmo 7 multietapa
    
    Gcdzs_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z); 

%% 31) Computar Normas H_2 e H_infinito para sistema projetado via 
%   Algoritmo 7

    %31.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_parcial=norm(Gcdys_parcial,2);
    
    %31.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norma_Hinf_Gcdzs_parcial=norm(Gcdzs_parcial,'inf');
    
%% 32) Computar Os polos do sistema via teorema 7

    Polos_Gcd_Parcial=cplxpair(pole(Gcdys_parcial));
    
%% 33) Obter Coeficientes de Amortecimento, Frequencia Natural e Medidas 
%   transitorias associadas

[wn_theo7,zeta_theo7,p_theo7]=damp(Gcdys_parcial);
cont=1;

%% 33) Computar medidas transitorias do sistema com retroalimentação projetada via algoritmo 5

Computar_Medidas_Transitorias_Malha_Fechada_via_parcial
    
%% 35) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_parcial

%% 36) Computar respostas de Gcdy(s) ao disturbio degrau unitario

    %36.1) Sinal y(t) gerado.
    
    [y_Gcdys_step_parcial,t_Gcdys_step_parcial,x_Gcdys_step_parcial]=...
        step(Gcdys_parcial,0:0.01:14);
    
    %36.2) Sinal u(t) gerado.
    
    u_Gcdys_step_parcial(1:size(x_Gcdys_step_parcial,1),1:n_u,1)=...
        x_Gcdys_step_parcial(:,:,1)*Kpf_parcial';
    
%% 35) Computar Tabelas

    Computar_Tabelas
          
%% 40) Construir graficos de Bode para sistema LCTI-MIMO em malha aberta e 
%   fechada.

    figure 
    subplot(133)
    sigmaplot(Gcdzs_classico,Gcdzs_parcial,Canon_Godzs)
    title('')
    legend('Gcdz(s) Classico','Gcdz(s) Parcial','Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    subplot(132)
    sigmaplot(Gcdzs_classico,Gcdzs_parcial)
    title('')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines
    
    subplot(131)
    sigmaplot(Gcdzs_classico)
    title('')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    sgtitle('Valores Singulares G_d_z(s)')
%% 41) Construir graficosde resposta transitoria de g(t)^2 ao disturbio  
%   impulso unitario do sistema LCTI-MIMO em malha aberta e fechada.

    figure
    subplot(133)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico,...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial,...
         t_Canon_Godys_impulse,y2_Canon_Godys_impulse);
    xlabel('t')
    ylabel('y^2(t)')
    legend('Gcdy(s) classico','Gcdy(s) parcial','Gody(s)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(132)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico,...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial);
    xlabel('t')
    ylabel('y^2(t)')
    grid on;
    xlim([0,7])
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(131)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico);
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0,2.5])
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Resposta y(t) ao impulso')
    
%% 42) Construir grafico de u(t) sistema LCTI-MIMO em malha aberta e 
%   fechada em resposta a disturbio degrau unitario.     
    
    figure
    subplot(121)
    plot(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,1),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,1));
    xlim([0 2.5])
    title('u^2_1(t)')
    xlabel('t')
    ylabel('u^2_1(t)')
    legend('G_c_d_y(s) classico','G_c_d_y(s) parcial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    plot(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,2),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,2));
    xlim([0 2.5])
    title('u^2_2(t)')
    xlabel('t')
    ylabel('u_2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('u(t) em resposta ao impulso')

%% 43) Graficos da resposta ao disturbio impulso da fracoes parciais e 
%      funcoes de transferencia.

    figure
    subplot(131)
    impulseplot(GAs,'-k',GNAs,'-r',Godys,'-b',Godzs,'-y');
    title('')
    legend('GA(s)','GNA(s)','Gody(s)','Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
    
    subplot(132)
    impulseplot(GAs,'-k',GNAs,'-r',Godys,'-b',Godzs,'-y');
    title('')
    grid on;
    xlim([0 75])
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
    
    subplot(133)
    impulseplot(GAs,'-k',GNAs,'-r',Godys,'-b');
    title('')
    grid on;
    xlim([0 75])
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
    
    sgtitle('Resposta ao impulso')
    
    figure
    sigma(GAs,'-k',GNAs,'-r',Godys,'-b',Godzs,'-y');
    title('Valores Singulares','FontSize',12)
    grid on;
    legend('GA(s)','GNA(s)','Gody(s)','Godz(s)','FontSize',10)
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
    
%% 44) Grafico para Localizacao dos Polos

    %44.1) Autovalores (Polos) a serem inseridos no grafico
    autovalores=[ Polos_Malha_Aberta',...
                  Polos_Gcd_classico',...  
                  Polos_Gcd_Parcial'];

    %44.2) Computando curvas de fronteira das regioes    
    [  x_horizontal_strip,y_inferior_horizontal_strip,...
        y_superior_horizontal_strip,...
        x_sector,y_inferior_sector,y_superior_sector,...
        x_parabola,y_inferior_parabola,y_superior_parabola,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = linhas_D_regioes(...
        autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

    %44.3) Grafico            
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
   ylim([-2. 2])
   xlim([-2.9 0])
   grid on
   sgtitle('Polos e localizacao Plano Complexo','FontSize',12)
