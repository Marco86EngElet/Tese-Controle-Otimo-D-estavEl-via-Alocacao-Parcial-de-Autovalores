clc, close, clear

format SHORT

%% 1) Definir quantidade de entradas e saidas.

    % 1.1) Quantidade de estados
    A_x=toeplitz([0.1:0.01:0.5,0.01]);
    n_x=size(A_x,1);
    
    % 1.2) Quantidade de saidas nao sem influencia direta de disturbios
    n_y=1;
    
    % 1.3) Quantidade de saidas com influencia direta de disturbios
    n_z=1;

    % 1.) Quantidade de de disturbios
    n_d=5;    

    % 1.) Quantidade de acoes de controle
    n_u=2;

    % ) Quantidade autovalores a serem modificados.
    n_ola=n_u;
 
%% 2) Inserir a quantidade de polos reais $n_r$ e par de polos complexos 
%   conjugados $n_{pcc}$.

    %2.1) Quantidade de polos reais
    n_r=2;
    
    %2.2) Quantidade de pares de polos complexos
    n_pcc=(n_x-2)/2;

%% 3) Delimitar limites para medidas transitoria associadas:

    %3.1) Sobressinal maximo.
    Mosmax=0.05;
    
    %3.2) Tempo de subida minimo e maximo.
    
    Trmin=0.74;
    Trmax=[];
    
    %3.3) Tempo de pico minimo.
    
    Tpmin=1.257;
    
    %3.4) Tempo de acomodacao maximo e minimo.
    
    csmax=0.01;
    Tsmax=4;
    Tsmin=2.4;
    
    %3.5) Tempo de atraso minimo.
    
    Tdmin=[];

%% 4) Computar os parametros das LMIs de D-estabilidade:

    Computar_parametros_D_estabilidade

%% 5) Representacao canonica das funcoes de transferencia do sistema em malha
%    Aberta
    
    %5.1) Criar matriz B_d
    B_d=0.1*eye(n_x,n_d);
    
    %5.2) Criar matriz C_y
    C_y=zeros(n_y,n_x);
    C_y(1,end)=1;
    C_y(1,end-1)=-0.5;
    C_y(1,end-2)=0.25;
    
    %5.3) Criar matriz C_z
    n_z=1;
    C_z=zeros(n_z,n_x);
    C_z(1,1)=1;
    C_z(1,end)=-1;
    
    %5.4) Criar matriz E_y
    E_y=zeros(n_y,n_d);
        
    %5.4) Criar matriz E_z
    E_z=zeros(n_z,n_d);
    E_z(1,1)=0.01;
    E_z(1,2)=-0.01;
    E_z(1,3)=0.01;
    E_z(1,4)=-0.01;
    E_z(1,5)=0.01;
    
    %5.5) Computa representacao canonica controlavel de Gody(s) 
    
    Canon_Godys = ss(A_x,B_d,C_y,E_y);
    
    %5.6) Computa representacao canonica controlavel de Godz(s)
    
    Canon_Godzs = ss(A_x,B_d,C_z,E_z);

%% 6) Computar Normas e Polos do sistema em Malha Aberta
    
    Computar_Normas_e_Polos_Malha_Aberta    
    
%% 7) Computar Frequencia Natural e Coeficiente de amortecimento

    [wn_God,zeta_God,p_God]=damp(Canon_Godys);
    
%% 8) Computar Medidas transitorias do sistema em Malha Aberta

    Computar_Medidas_Transitoria_Sistema_Malha_Aberta

%% 9) Computar respostas de Gody(s) ao disturbio impulso unitario

    tempo_simulacao_impulso=0:0.03:30;
    
    Computar_Dados_Simulacao_impulso
    
%% 10) Especificar as Matrizes associadas ao controle 'B_u', 'D_y', 'D_z'.
    
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
    
    c_H2=[];
    c_Hinf=[];

%% 13) Executar Algoritmo 2 

    %13.1) Iniciar matriz de retroalimentacao 
    
    Kpf_parcial=zeros(n_u,n_x);

    %13.2) Entrada Q_p
    
    Q_p = sqrt(2)*eye(2);

    %13.3) Obter Matrizes dos Autovalores e Autovetores
        
    [Right_Eigenvectors, Eigenvalues, Left_Eigenvectors] = ...
        eig(A_x + B_u * Kpf_parcial);
        
    %13.4) Encontrar as matrizes para alocacao
    [row_poles,col_poles]=find(real(Eigenvalues)>0);
    
    Lambda_p = Eigenvalues(row_poles, col_poles);
    L_p = Left_Eigenvectors(1:end, col_poles);    
    
    %13.4) Computar controle otimo misto D-estavel com alocacao parcial
    Controle_otimo_misto_D_estavel_alocacao_parcial

        %13.5) Iniciar contagem de tempo para otimizacao multietapa via 
        %       Algoritmo 7
        tic
    
        %13.6) Executar Algoritmo 7 Multietapa   
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
    
        %13.7) Extrair variaveis de decisao para computar K_pf 

        tempo_otimizacao_classico=toc;
        
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X);    

%% 14) Computar Matriz de retroalimentacao via Algoritmo 3

Kpf_classico=optimal_Tilde_W/optimal_Tilde_X;
Kpf_classico=Kpf_classico* (Q_p * L_p')/2;

%% 15) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 3

    % 15.1) Espaco de Estados para Gcdy(s) 
    
    Gcdys_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_y+D_y*Kpf_classico,E_y);

    % 15.2) Espaco de Estados para Gcdz(s)  
    
    Gcdzs_classico=...
        ss(A_x+B_u*Kpf_classico,B_d,C_z+D_z*Kpf_classico,E_z);

%% 16) Computar medidas transitorias e tabela para sistema em malha fechada projetado via algoritmo 3 

    Computar_Medidas_Transitorias_Malha_Fechada_via_classico
    
%% 17) Computar Normas H_2 e H_infinito para sistema projetado via 
%   Algoritmo 3

    %17.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_classico=norm(Gcdys_classico,2);
    
    %17.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norm_Hinf_Gcdzs_classico=norm(Gcdzs_classico,'inf');
    
%% 18) Computar Os polos do sistema projetado via Algoritmo 3

    Polos_Gcd_classico=cplxpair(pole(Gcdys_classico));
    
%% 19) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_livre

%% 20) Iniciar contagem de tempo para otimizacao multietapa via 
%         Algoritmo 7
        
        c_H2=1;
        c_Hinf=2;
        tic
    
        %20.1) Executar Algoritmo 7 Multietapa   
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
        
        optimal_Tilde_W = value(Tilde_W);
        optimal_Tilde_X = value(Tilde_X); 
        
        tempo_otimizacao_parcial = toc;
        
        KD_parcial=optimal_Tilde_W*inv(optimal_Tilde_X);
        Kpf_parcial = Kpf_parcial + KD_parcial * (Q_p * L_p') / 2;

    % 20.2) Finalizar tempo de otimizacao multietapa
    
%% 21) Representacao em espaco de estados do sistema em malha fechada 
%   projetado via Algoritmo 7 multietapa

    %21.1) Espaco de estados para Gcdy(s) via algoritmo 7 multietapa
    
    Gcdys_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);

    %21.2) Espaco de estados para Gcdz(s) via algoritmo 7 multietapa
    
    Gcdzs_parcial=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z); 

%% 22) Computar Normas H_2 e H_infinito para sistema projetado via 
%   Algoritmo 7

    %22.1) Norma H_2 do sistema canonico controlavel de Gcdy(s)
    
    norma_H2_Gcdys_parcial=norm(Gcdys_parcial,2);
    
    %22.2) Norma H_infinito do sistema canonico controlavel de Gcdz(s)
    
    norma_Hinf_Gcdzs_parcial=norm(Gcdzs_parcial,'inf');
    
%% 23) Computar Os polos do sistema via teorema 7

    Polos_Gcd_Parcial=cplxpair(pole(Gcdys_parcial));
    
%% 24) Obter Coeficientes de Amortecimento, Frequencia Natural e Medidas 
%   transitorias associadas

[wn_theo7,zeta_theo7,p_theo7]=damp(Gcdys_parcial);
cont=1;

%% 25) Computar medidas transitorias do sistema com retroalimentação projetada via algoritmo 5

Computar_Medidas_Transitorias_Malha_Fechada_via_parcial
    
%% 26) Computar respostas de Gcdy(s) ao disturbio impulso unitario

    Computar_yt_ut_resposta_ao_impulso_metodo_parcial

%% 27) Computar Tabelas

    Computar_Tabelas
          
%% 28) Construir graficos de Bode para sistema LCTI-MIMO em malha aberta e 
%   fechada.

    figure 
    subplot(121)
    sigmaplot(Gcdzs_classico,Gcdzs_parcial,Canon_Godzs)
    title('')
    legend('Gcdz(s) Classico','Gcdz(s) Parcial','Godz(s)','FontSize',10)
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    subplot(122)
    sigmaplot(Gcdzs_classico,Gcdzs_parcial)
    title('')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); % Set thickness for all lines

    sgtitle('Valores Singulares G_d_z(s)')
    
%% 29) Construir graficosde resposta transitoria de g(t)^2 ao disturbio  
%   impulso unitario do sistema LCTI-MIMO em malha aberta e fechada.

    figure
    
    subplot(151)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico(1:end,1,1),...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial(1:end,1,1));
    title('y^2(t) X d_1(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(152)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico(1:end,1,2),...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial(1:end,1,2));
    title('y^2(t) X d_2(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(153)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico(1:end,1,3),...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial(1:end,1,3));
    title('y^2(t) X d_3(t)')
    xlim([0 7])
    xlabel('t')
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(154)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico(1:end,1,4),...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial(1:end,1,4));
    title('y^2(t) X d_4(t)')
    xlabel('t')
    xlim([0 7])
    ylabel('y^2(t)')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(155)
    plot(t_Gcdys_impulse_classico,y2_Gcdys_impulse_classico(1:end,1,5),...
         t_Gcdys_impulse_parcial,y2_Gcdys_impulse_parcial(1:end,1,5));
    title('y^2(t) X d_5(t)')
    xlabel('t')
    ylabel('y^2(t)')
    xlim([0 7])
    legend('Gcdy(s) classico','Gcdy(s) parcial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Resposta y^2(t) ao impulso')
    
%% 30) Construir grafico de u(t) sistema LCTI-MIMO em malha aberta e 
%   fechada em resposta a disturbio degrau unitario.     
    
    figure
    subplot(121)
    semilogy(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,1),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,1));
    title('u_1(t) em resposta ao degrau')
    xlabel('t')
    ylabel('u^2_1(t)')
    xlim([0,20]);
    legend('G_c_d_y(s) classico','G_c_d_y(s) parcial')
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    subplot(122)
    semilogy(t_Gcdys_impulse_classico,u2_Gcdys_impulse_classico(:,2),...
         t_Gcdys_impulse_parcial,u2_Gcdys_impulse_parcial(:,2));
    title('u^2_2(t) em resposta ao degrau')
    xlabel('t')
    ylabel('u_2(t)')
    xlim([0,20]);
    grid on;
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3); 
    
    sgtitle('Resposta u(t) ao degrau')

%% 31) Grafico para Localizacao dos Polos

    %31.1) Autovalores (Polos) a serem inseridos no grafico
    autovalores=[ Polos_Malha_Aberta',...
                  Polos_Gcd_classico',...  
                  Polos_Gcd_Parcial'];

    %31.2) Computando curvas de fronteira das regioes    
    [  x_horizontal_strip,y_inferior_horizontal_strip,...
        y_superior_horizontal_strip,...
        x_sector,y_inferior_sector,y_superior_sector,...
        x_parabola,y_inferior_parabola,y_superior_parabola,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = linhas_D_regioes(...
        autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

    %31.3) Grafico            
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
   ylim([-1.1 1.1])
   xlim([-2 0])
   grid on
   sgtitle('Polos e localizacao Plano Complexo','FontSize',12)
   
