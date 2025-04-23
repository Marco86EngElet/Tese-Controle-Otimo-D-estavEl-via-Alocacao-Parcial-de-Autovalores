clc,close all,clear all,
%   1) Criar matriz A_x formato de toeplitz
    A_x=toeplitz([0.1:0.01:0.5,0.01]);
    n_x=size(A_x,1);
    %---------------------------------
%   2)  Criar matriz B_u
    n_u=2;
    B_u=zeros(n_x,n_u);
    B_u(end-1:end,1:2)=eye(2);
    %---------------------------------
%   3) Criar matriz B_d
    n_d=5;
    B_d=0.1*eye(n_x,n_d);
    %---------------------------------
%   4) Criar matriz C_y
    n_y=1;
    C_y=zeros(n_y,n_x);
    C_y(1,end)=1;
    C_y(1,end-1)=-0.5;
    C_y(1,end-2)=0.25;
    %---------------------------------
%   5) Criar matriz D_y
    D_y=zeros(n_y,n_u);
    D_y(1,1)=0.01;
    %--------------------------------    
%   6) Criar matriz E_y
    E_y=zeros(n_y,n_d);
    %--------------------------------
%   7) Criar matriz C_z
    n_z=1;
    C_z=zeros(n_z,n_x);
    C_z(1,1)=1;
    C_z(1,end)=-1;
    %-----------------------------  
%   8) Criar matriz D_z
    D_z=zeros(n_z,n_u);
    D_z(1,1)=0.01;
    %----------------------------
%   9) Criar matriz E_z
    E_z=zeros(n_z,n_d);
    E_z(1,1)=0.01;
    E_z(1,2)=-0.01;
    E_z(1,3)=0.01;
    E_z(1,4)=-0.01;
    E_z(1,5)=0.01;
    %------------------------------
%   10) Computar quantidade de polos instaveis
    quantidade_polos_instaveis=sum(eig(A_x)>0);
    n_ola=quantidade_polos_instaveis;
    %-----------------------
%   11) Localizar autovalores instaveis
    [n_l,n_c]=find(eig(A_x)>0);
    %-------------------------
%   12)  Computar matrizes associadas aos polos de Gd e necessarias para o 
%       metodo de alocacao
    [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=eig(A_x);
    Lambda= ...
        Eigenvalues(n_l,n_l);
    L_j = Left_Eigenvectors(1:end,n_l);
    Q = sqrt(2)*eye(n_ola);      
    %-----------------------
%   12) Matriz de controlabilidade parcial (Sistema Reduzido)
    matriz_controlabilidade_parcial=ctrb(Lambda,(Q*L_j')*B_u/2);
    %---------------------
%   13) Quantidade de estados nao controlaveis para sistema reduzido
    estados_nao_controlaveis_sistema_reduzido=...
        n_ola-rank(matriz_controlabilidade_parcial);   
    %-----------------------
%   14) Caracteristicas da regiao de D-estabilidade
        %--------------------
        %14.1) Angulo dos polos entre '-theta_s' e '+theta_s'
        mosmax=0.05;
        theta_s=atan(pi/log(mosmax^(-1)));
        %------------------
        %14.2) Parte real inferior a -alpha_v
        Tsmax=4;
        csmax=0.01;
        av1=[];
        av2=4/Tsmax;
        av3=(log(csc(theta_s)*csmax^(-1)))/Tsmax;
        alpha_v=max([av1,av2,av3]);
        %-------------------
        %14.3) Parte real superior a -beta_v
        Trmin=0.74;
        Tsmin=2.4;
        bv1=2.2/Trmin;
        bv2=4/Tsmin;
        beta_v=min([bv1,bv2]);
        %-------------------
        %14.4) Raio disco
        r_d=(1.1+0.125*cos(theta_s)+0.495*cos(theta_s)^2)/Trmin;;
        %------------------
        %14.5) Centro Disco
        q_d=0;
        %-----------------
        %14.6) Parte Imaginaria entre '-w_H' e "+w_H"
        Tpmin=1.257;
        w_H=pi/Tpmin;
        %---------------
        %14.7) Fator de amortecimento Parabola de Estabilidade
        e_P=[]; 
%----------------------
%   15) Configuracoes para YALMIP
        %15.1) Configuracoes para otimizacao semidefinida no YALMIP
        Yalmip_sdpsettings =...
             sdpsettings('verbose',1,'solver','lmilab','debug',1);
        %-------------
        %15.2) Pesos da funcao custo de otimizacao 
        c_H2=1;
        c_Hinf=2;
        %-------------
%   16) Computar matriz de retroalimentacao com alocacao parcial
    if estados_nao_controlaveis_sistema_reduzido==0
        %16.1) Obter matrizes do sistema reduzido
        [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz]=...
           matrizes_sistema_reduzido(Lambda,L_j,Q,B_u,B_d,C_y,C_z);
        %-------------------
        %16.2) Executar metodo de otimizacao para alocacao parcial otima
        tic
        controle_misto_H2_Hinf_Parcial_D_estabilidade
        tempo_alocacao_otima=toc;
        %---------------------
        %16.3) Computar matriz de retroalimentacao otima
        Kpf_parcial=KD_parcial*(Q*L_j')/2;
        %------------------
        %16.4) Executar metodo de otimizacao para alocacao parcial factivel
        %   D-estavel
        C_H2=[];
        c_Hinf=[];
        tic
        controle_misto_H2_Hinf_Parcial_D_estabilidade
        tempo_alocacao_destavel=toc;
        Kpf_destavel=KD_parcial*(Q*L_j')/2;
    else
       disp('sistema reduzido nao e controlavel') 
    end      
    %-----------------------------
%=================================        
%   17) Objeto para sistema em malha aberta para Gdy e Gdz
    sistema_malha_aberta_Gdy=ss(A_x,B_d,C_y,E_y);
    sistema_malha_aberta_Gdz=ss(A_x,B_d,C_z,E_z);    
    %-------------------------
%   18) Objeto para sistema em malha fechada para Gdy e Gdz com otimizacao
    %18.1) Sistema d(t) para y(t)  
    sistema_alocacao_parcial_Gdy=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);
    %----------------------------
    %18.2) Sistema d(t) para z(t) 
    sistema_alocacao_parcial_Gdz=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z); 
    %----------------------------
%   19) Objeto para sistema em malha fechada para Gdy e Gdz so com
%           D-estabilidade
    %19.1) Sistema d(t) para y(t)  
    sistema_alocacao_destavel_Gdy=...
        ss(A_x+B_u*Kpf_destavel,B_d,C_y+D_y*Kpf_destavel,E_y);
    %---------------------------------    
    %19.2) Sistema d(t) para z(t)  
    sistema_alocacao_destavel_Gdz=...
        ss(A_x+B_u*Kpf_destavel,B_d,C_z+D_z*Kpf_destavel,E_z); 
    %--------------------------------    
%   20) Tabela dos Polos  
    Polos_Malha_Aberta=...
        cplxpair(pole(sistema_malha_aberta_Gdy));

    Polos_Malha_Fechada_Alocacao_Parcial_Otima=...
        cplxpair(pole(sistema_alocacao_parcial_Gdy));
    
    Polos_Malha_Fechada_Alocacao_Parcial_Destavel=...
        cplxpair(pole(sistema_alocacao_destavel_Gdy));
    
    Tabela_Polos=...
        table( Polos_Malha_Aberta,...
               Polos_Malha_Fechada_Alocacao_Parcial_Otima);
    %----------------------------    
%   21) Graficos para Polos e posicao

    autovalores=...
        [ Polos_Malha_Aberta',...
          Polos_Malha_Fechada_Alocacao_Parcial_Otima',...
          Polos_Malha_Fechada_Alocacao_Parcial_Destavel'];
          
    [  x_horizontal_strip,y_inferior_horizontal_strip,...
            y_superior_horizontal_strip,...
            x_sector,y_inferior_sector,y_superior_sector,...
            x_parabola,y_inferior_parabola,y_superior_parabola,...
            x_disk,y_disk,...
            x_alphav,y_alphav,...
            x_betav,y_betav ] = linhas_D_regioes(...
            autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P);  
        
    figure
    subplot(121)
    plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
          real(Polos_Malha_Fechada_Alocacao_Parcial_Destavel),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial_Destavel),...
          '+r',real(Polos_Malha_Fechada_Alocacao_Parcial_Otima),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial_Otima),'xb',...
          x_alphav,y_alphav,'-k',...
          x_betav,y_betav,'-k',...
          x_disk,y_disk,'-k',...
          x_horizontal_strip,y_inferior_horizontal_strip,'-k',...
          x_horizontal_strip,y_superior_horizontal_strip,'-k',...
          x_sector,y_inferior_sector,'-k',...
          x_sector,y_superior_sector,'-k',...
          x_parabola,y_inferior_parabola,'-k',...
          x_parabola,y_superior_parabola,'-k',...
          'LineWidth',2,'MarkerSize',8)
          xlim([-2 11]);
 title('Polos e localizacao Plano Complexo','FontSize',12)
 legend('Aberta','Destavel','Otima','Fronteira','FontSize',12)
 
 subplot(122)
    plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
          real(Polos_Malha_Fechada_Alocacao_Parcial_Destavel),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial_Destavel),'+r',...
          real(Polos_Malha_Fechada_Alocacao_Parcial_Otima),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial_Otima),'xb',...
          x_alphav,y_alphav,'-k',...
          x_betav,y_betav,'-k',...
          x_disk,y_disk,'-k',...
          x_horizontal_strip,y_inferior_horizontal_strip,'-k',...
          x_horizontal_strip,y_superior_horizontal_strip,'-k',...
          x_sector,y_inferior_sector,'-k',...
          x_sector,y_superior_sector,'-k',...
          x_parabola,y_inferior_parabola,'-k',...
          x_parabola,y_superior_parabola,'-k',...
          'LineWidth',2,'MarkerSize',8)
          xlim([-2 0.1]);

 
 %-------------------------------
 %  22) Tabela para comparar normas dos sitemas

    %22.1) Computar normas do sistema malha aberta

    norma_H2(1,1)=norm(sistema_malha_aberta_Gdy,2);
    
    variacao_percentual_norma_H2(1,1)=0;
    
    norma_Hinf(1,1)=norm(sistema_malha_aberta_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(1,1)=0;
    %----------------------------        
    %22.2) Computar normas do sistema malha fechada alocacao livre

    norma_H2(2,1)=norm(sistema_alocacao_destavel_Gdy,2);
    
    variacao_percentual_norma_H2(2,1)=...
        100*(norma_H2(2,1)-norma_H2(1,1))/norma_H2(1,1);
    
    norma_Hinf(2,1)=norm(sistema_alocacao_destavel_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(2,1)=...
        100*(norma_Hinf(2,1)-norma_Hinf(1,1))/norma_Hinf(1,1);
    %----------------------------
        
    %22.3) Computar normas do sistema malha fechada alocacao parcial

    norma_H2(3,1)=norm(sistema_alocacao_parcial_Gdy,2);
    
    variacao_percentual_norma_H2(3,1)=...
        100*(norma_H2(3,1)-norma_H2(1,1))/norma_H2(1,1);
    
    norma_Hinf(3,1)=norm(sistema_alocacao_parcial_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(3,1)=...
        100*(norma_Hinf(3,1)-norma_Hinf(1,1))/norma_Hinf(1,1);
    %----------------------------        
    %22.4) Criar tabela

    Funcoes_transferencia=["Aberta";"Factivel Destavel";...
                            "Alocacao Otima"];
    Tabela_Normas = table(Funcoes_transferencia,...
        norma_H2,variacao_percentual_norma_H2,...
        norma_Hinf,variacao_percentual_norma_Hinf),
%================================
%   23) Graficos
%   23.1) 
 time_opt=timeoptions;
    time_opt.InputLabels.FontSize=12;
    time_opt.OutputLabels.FontSize=12;
    time_opt.XLabel.FontSize=12;
    time_opt.YLabel.FontSize=12;
    time_opt.TickLabel.FontSize=12;

    bode_options=bodeoptions;
    bode_options.PhaseVisible='off';
    bode_options.InputLabels.FontSize=12;
    bode_options.OutputLabels.FontSize=12;
    bode_options.Title.FontSize=12;
    bode_options.XLabel.FontSize=12;
    bode_options.YLabel.FontSize=12;
    bode_options.TickLabel.FontSize=12;
    %--------------------------------
%   23.2) Resposta degrau y(t)        
    figure
        stepplot(sistema_alocacao_destavel_Gdy,'-r',...
             sistema_alocacao_parcial_Gdy,'-b',time_opt)
        title('Resposta y(t) ao degrau','FontSize',12)
        legend('Factivel D-estavel','Parcial otimo','FontSize',12)

        % Encontrar todos os eixos da figura
        axesArray = findobj(gcf, 'Type', 'Axes');  
        for k = 1:length(axesArray)
            xlim(axesArray(k), [0 900]);  % Aplica xlim a cada eixo
        end
        %----------------------------

    figure    
        stepplot(sistema_alocacao_destavel_Gdz,'-r',...
             sistema_alocacao_parcial_Gdz,'-b',time_opt)
        title('Resposta z(t) ao degrau','FontSize',12)
        legend('Factivel D-estavel','Parcial otimo','FontSize',12)

        % Encontrar todos os eixos da figura
        axesArray = findobj(gcf, 'Type', 'Axes');  
        for k = 1:length(axesArray)
            xlim(axesArray(k), [0 900]);  % Aplica xlim a cada eixo
        end
        %---------------------

    figure
        bodemag(sistema_alocacao_destavel_Gdy,'-r',...
             sistema_alocacao_parcial_Gdy,'-b',bode_options)
        title('Gdy(jw)','FontSize',12)
        legend('Factivel D-estavel','Parcial otimo','FontSize',12)
        %------------------

    figure    
        bodemag(sistema_alocacao_destavel_Gdz,'-r',...
             sistema_alocacao_parcial_Gdz,'-b',bode_options)
        title('Gdz(jw)','FontSize',12)
        legend('Factivel D-estavel','Parcial otimo','FontSize',12)
        %----------------------

  