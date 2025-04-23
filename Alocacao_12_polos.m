clear all,close all,clc,
%1)  Criando Funcoes de Transferencia Segunda Ordem 
%    Polos nao dominante
%
%   1.1) Especifica-se o criterio de acomodao 2%   
    csmax=0.02;
%   --------------------
%   1.2) Inicializar funcao GNA(s)
    GNA=0;
%   ---------------------
%   1.3) Inicializar funcao Gody(s)
    Gody=0;
%   ---------------------
%   1.4) Inicializar funcao Godz(s)
    Godz=1;
%   -----------------------
%   1.5) Quantidade de polos nao alocaveis
    num_pna=5;
%   -----------------------
%   1.6) Iniciar contagem para computar fracoes parciais de GNA(s)
for i=1:num_pna
    %1.6.1) coeficiente de amortecimento fracao parcial G(i,1)
    zeta(i,1)=exp(-2+i/30);
    %----------------------
    %1.6.2) frequencia natural fracao parcial G(i,1)
    wn(i,1)=0.9^(i/30);
    %----------------------
    %1.6.3) parte real do polo da fracao parcial G(i,1)
    parte_real(i,1)=wn(i,1)*zeta(i,1);
    %-----------------------
    %1.6.4) parte imaginaria do polo da fracao parcial G(i,1)
    parte_imaginaria(i,1)=wn(i,1)*(sqrt(1-zeta(i,1)^2));
    %-----------------------
    %1.6.5) ganho do polo da fracao parcial G(i,1)
    ganho(i,1)=0.15;
    %-----------------------
    %1.6.6) fracao parcial G(i,1)
    G(i,1)=...
        Criar_Funcao_Transferencia_SISO_2nd(ganho(i,1),...
            zeta(i,1),wn(i,1));
        %-----------------------
    %1.6.7) Medidas transitorias associadas a G(i,1)
    [mos(i,1),ts(i,1),td(i,1),tr(i,1),tp(i,1)] =...
        medidas_SISO_2nd(zeta(i,1),wn(i,1),csmax);
    %------------------------------
    %1.6.8) Soma das fracoes parciais nao alocaveis
    GNA=GNA+G(i,1);    
    %-------------------------
    %1.6.9) adicionar G(i,1) a funcao de transferencia d(t) para y(t)
    Gody=Gody+G(i,1);
    %--------------------------
    %1.6.10) adicionar G(i,1) a funcao de transferencia d(t) para y(t)
    Godz=Godz+G(i,1);
    %-------------------------
end
%   1.7) Computar maximo valor de cada G(i,1)
    max_value=(1+mos).*ganho;
    %--------------------------
%   1.8) Tabela para obter medidas transitorias de G(i,1) nao alocaveis
    Tabela_Medidas_Transitorias_Gi_nao_Alocaveis=...
        table(zeta,wn,parte_real,parte_imaginaria,...
        max_value,mos,ts,td,tr,tp);
    %------------------------------------
%==============================
%   2) Criando Funcao GA
%   2.1) frequencia natural
    wn(num_pna+1,1)=0.2;
    %--------------------------
%   2.2) coeficiente de amortecimento
    zeta(num_pna+1,1)=0.1;
    %-------------------------
%   2.3) ganho
    ganho(num_pna+1,1)=7;
    %--------------------------
%   2.4) parte real do polo
    parte_real(num_pna+1,1)=...
        wn(num_pna+1,1)*zeta(num_pna+1,1);
    %---------------------------
%   2.5) parte imaginaria do polo
    parte_imaginaria(num_pna+1,1)=...
        wn(num_pna+1,1)*(sqrt(1-zeta(num_pna+1,1)^2));
    %-----------------------------
%   2.6) criar funcao Gd (Funcao de transferencia com dinamica dominante)
    GA=Criar_Funcao_Transferencia_SISO_2nd(...
        ganho(num_pna+1,1),zeta(num_pna+1,1),wn(num_pna+1,1));
%==============================
%   3) Computar funcao Gody(s) e Godz(s)
%   3.1) Medidas transitorias associada a Gd(s)
    [mos(num_pna+1,1),ts(num_pna+1,1),...
        td(num_pna+1,1),tr(num_pna+1,1),tp(num_pna+1,1)] =...
        medidas_SISO_2nd(zeta(num_pna+1,1),wn(num_pna+1,1),csmax);
    %-----------------------
%   3.2) Adicionar maximo valor de Gd(s) a 'max_value'
    max_value=mos.*ganho;
    %----------------------
%   3.3) adicionar GA(s) a funcao Gody(s) 
    Gody=Gody+GA;
    %--------------------
%   3.4) adicionar GA(s) a funcao Godz(s)
    Godz=Godz+GA;
    %---------------------
%=========================
%   4) Obter tabela para medidas transitorias de todas as funcoes parciais
    Tabela_Funcoes_Parciais=...
        table(zeta,wn,parte_real,parte_imaginaria,...
            max_value,mos,ts,td,tr,tp);
%=========================
%   5) Verificar resposta ao degrau e graficos de bode para Gd(s), Gnd(s),
%   Gdy(s) e Gdz(s)
%   5.1) Configuracoes para grafico de resposta a entrada degrau d(t)
    time_opt=timeoptions;
    time_opt.InputLabels.FontSize=12;
    time_opt.OutputLabels.FontSize=12;
    time_opt.XLabel.FontSize=12;
    time_opt.YLabel.FontSize=12;
    time_opt.TickLabel.FontSize=12;
    %---------------------
%   5.2) Configuracoes para grafico de bode    
    bode_options=bodeoptions;
    bode_options.PhaseVisible='off';
    bode_options.InputLabels.FontSize=12;
    bode_options.OutputLabels.FontSize=12;
    bode_options.Title.FontSize=12;
    bode_options.XLabel.FontSize=12;
    bode_options.YLabel.FontSize=12;
    bode_options.TickLabel.FontSize=12;
    %-------------------------
%   5.3) Graficos    
    figure
    subplot(121)
    stepplot(GA,'-r',GNA,'-y',Gody,'-b',Godz,'-m',time_opt);
    legend('GA(s)','GNA(s)','Gody(s)','Godz(s)','FontSize',12)
    title('Resposta degrau','FontSize',12)
    subplot(122)
    bodeplot(GA,'-r',GNA,'-y',Gody,'-b',Godz,'-m',bode_options);
    legend('GA(s)','GNA(s)','Gody(s)','Godz(s)','FontSize',12)
    title('Bodemag','FontSize',12)
%===========================
%   6) Computar modelos de espaco de estados em malha aberta
%   6.1) Modelo Espaco de estados de Gody(s)
    sistema_malha_aberta_Gdy = ss(Gody);
    %---------------------
%   6.2) Modelo Espaco de estados de Gxz
    sistema_malha_aberta_Gdz = ss(Godz);
    %------------------------
%   6.3) Computar as matrizes do modelo espaco de estados
    A_x=sistema_malha_aberta_Gdy.A;
    B_d=sistema_malha_aberta_Gdy.B;
    C_y=sistema_malha_aberta_Gdy.C;
    C_z=sistema_malha_aberta_Gdz.C;
    E_y=sistema_malha_aberta_Gdy.D;
    E_z=sistema_malha_aberta_Gdz.D;    
    %-----------------------
%   6.4) Obter dimensoes dos vetores de estados e saidas
    n_x=size(A_x,1);
    n_y=size(C_y,1);
    n_z=size(C_z,1);
    n_d=size(B_d,2);
    %----------------------
%==============================
%   7) Computar matrizes para alocacao parcial
%   7.1) Obter Autovalores e Autovetores de A_x
    [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=eig(A_x);
    %---------------------------
%   7.2) Encontrar posicao autovalores e autovetores associados ao polos de
%       Gd(s)
    [posicao_linhas_polos_Gd,...
    posicao_colunas_polos_Gd]=...
        find(real(Eigenvalues)<-0.019 & real(Eigenvalues)>-0.021);
    %-------------------------
%   7.3) Computar matrizes associadas aos polos de Gd e necessarias para o 
%       metodo de alocacao
    Lambda=Eigenvalues(...
        posicao_linhas_polos_Gd,posicao_colunas_polos_Gd);
    L_j = Left_Eigenvectors(1:end,posicao_colunas_polos_Gd);
    Q = [1 1;1i -1i]; 
    %-----------------------
%   7.4) Obter dimensao do vetor das variaveis de controle para alocacao 
%       parcial
    n_ola=size(Lambda,1);
    n_u=n_ola;
    %-----------------------
%   7.5) Espeficar B_u, D_u2, D_uinf 
    B_u=zeros(n_x,n_u);
    B_u(end-1,1)=1;
    B_u(end,2)=-1;
    D_y=zeros(n_y,n_u);
    D_y(1,1)=0.1;
    D_z=D_y;
%========================
%   8) Checar quantidade de estados nao controlaveis dos modelos
%   8.1) Matriz de Controlabilidade Total
    matriz_controlabilidade_total=ctrb(A_x,B_u);
    %-----------------
%   8.2) Matriz de Controlabilidade Total     
    %Quantidade de estados nao controlaveis para sistema original
    estados_nao_controlaveis_sistema_original=...
        n_x-rank(matriz_controlabilidade_total);
    %-----------------
%   8.3) Matriz de controlabilidade parcial (Sistema Reduzido)
    matriz_controlabilidade_parcial=ctrb(Lambda,(Q*L_j')*B_u/2);
    %------------------
%   8.4) Quantidade de estados nao controlaveis para sistema reduzido
    estados_nao_controlaveis_sistema_reduzido=...
        n_ola-rank(matriz_controlabilidade_parcial);
%========================
%   9) Caracteristicas da regiao de D-estabilidade
%       ------------------
%       9.1) Parte Imaginaria entre '-w_H' e "+w_H"
        %Tpmin=3.15;
        %w_H=pi/Tpmin;
        w_H=[];
%       ----------------
%       9.2) Angulo dos polos entre '-theta_s' e '+theta_s'
        mosmax=0.06;
        theta_s=atan(pi/log(mosmax^(-1)));
%       ---------------
%       9.3) Raio disco
       Trmin=0.5;
        r_d=(1.1+0.125*cos(theta_s)+0.495*cos(theta_s)^2)/Trmin;
%        r_d=[];
%        ---------------
%       9.4) Parte real inferior a '-alpha_v'
%        Trmax=16;
        Tsmax=4.71;
%       av1=2.2/Trmax;
        av1=[];
        av2=4/Tsmax;
        csmax=0.02;
        av3=(log(csc(theta_s)*csmax^(-1)))/Tsmax;
        alpha_v=max([av1,av2,av3]);
%       -------------------
%       9.5) Parte real superior a '-beta_v'
%        Tsmin=6;
        bv1=2.2/Trmin;
%        bv1=[];
%        bv2=4/Tsmin;
        bv2=[];
        beta_v=min([bv1,bv2]);
%       ---------------------
%       9.6) Centro Disco
        q_d=0;
%        -------------------
%       9.7) Fator de amortecimento Parabola de Estabilidade
        e_P=[];
%       -------------------------
%============================
%   10)  Configuracoes para YALMIP
%   10.1) Configuracoes para otimizacao semidefinida no YALMIP
        Yalmip_sdpsettings =...
             sdpsettings('verbose',1,'solver','lmilab','debug',1);
    %------------------
%   10.2) Pesos da funcao custo de otimizacao 
        c_H2=2;
        c_Hinf=8;
%=================
%   11) Computar matriz de retroalimentacao com alocacao livre
    tic
    controle_misto_H2_Hinf_D_estabilidade
    tempo_otimizacao_alocacao_livre=toc;
%===============
%   12) Computar matriz de retroalimentacao com alocacao parcial
%   12.1) Construir sistema reduzido em malha aberta para alocacao parcial
    if estados_nao_controlaveis_sistema_reduzido==0
%       12.1.1) Obter matrizes do sistema reduzido
        [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz]=...
         matrizes_sistema_reduzido(Lambda,L_j,Q,B_u,B_d,C_y,C_z);
        %------------
%       12.1.2) Executar metodo de otimizacao para alocacao parcial
        tic
        controle_misto_H2_Hinf_Parcial_D_estabilidade
        tempo_otimizacao_alocacao_parcial=toc;
        %-------------
%       12.1.3) Computar matriz de retroalimentacao
        Kpf_parcial=KD_parcial*(Q*L_j')/2;
        %------------
    else
       disp('sistema reduzido nao e controlavel') 
    end
%===================
%   13) Sistema malha fechada via controle otimo alocacao livre
%   13.1) Sistema d(t) para y(t)  
    sistema_alocacao_livre_Gdy=...
        ss(A_x+B_u*Kpf_livre,B_d,C_y+D_y*Kpf_livre,E_y);
    %---------------
%   13.2) Sistema d(t) para z(t)  
    sistema_alocacao_livre_Gdz=...
        ss(A_x+B_u*Kpf_livre,B_d,C_z+D_z*Kpf_livre,E_z);
%===============
%   14) Sistema malha fechada via controle otimo alocacao parcial
%   14.1) Sistema d(t) para y(t)  
    sistema_alocacao_parcial_Gdy=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);
    %----------
%   14.2) Sistema d(t) para z(t)  
    sistema_alocacao_parcial_Gdz=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z);   
    %--------------
%==============
%   15) Graficos
%   15.1) Comparar resposta transitoria a disturbio grau unitario em y(t)
    figure
    subplot(121)
    stepplot(sistema_alocacao_livre_Gdy,'-r',...
         sistema_alocacao_parcial_Gdy,'-b',...
         sistema_malha_aberta_Gdy,'-k',time_opt)
    title('Resposta y(t) ao degrau','FontSize',12)
    legend('Livre','Parcial','Gody(s)','FontSize',12)
    subplot(122)
    stepplot(sistema_alocacao_parcial_Gdy,'-b',GNA,'-m',time_opt)
    title('Resposta y(t) ao degrau','FontSize',12)
    legend('Parcial','GNA(s)','FontSize',12)
    
    figure
    subplot(121)
    stepplot(sistema_alocacao_livre_Gdz,'-r',...
         sistema_alocacao_parcial_Gdz,'-b',...
         sistema_malha_aberta_Gdz,'-k')
    title('Resposta z(t) ao degrau')
    legend('Livre','Parcial','Godz(s)','FontSize',12)
    subplot(122)
    stepplot(sistema_alocacao_parcial_Gdz,'-b',GNA,'-m')
    title('Resposta z(t) ao degrau','FontSize',12)
    legend('Parcial','GNA(s)','FontSize',12)    
    %-----------------
%   15.2)  Comparar graficos de bode para Gxy(s)
    figure
    subplot(121)
    bodeplot(sistema_alocacao_livre_Gdy,'-r',...
            sistema_alocacao_parcial_Gdy,'-b',...
            sistema_malha_aberta_Gdy,'-k',...
            GNA,'-m',bode_options)
    legend('Livre','Parcial','Gody(s)','Aberto GNA(s)','FontSize',12)
    title('Bodemag Gdy(j\omega)','FontSize',12)
    subplot(122)
    bodeplot(sistema_alocacao_livre_Gdz,'-r',...
            sistema_alocacao_parcial_Gdz,'-b',...
            sistema_malha_aberta_Gdz,'-k',...
            GNA,'-m',bode_options)
    legend('Livre','Parcial','Godz(s)','GNA','FontSize',12)
    title('Bodemagz Gdz(j\omega)','FontSize',12)
    %--------------
%   15.3) Tabela dos Polos  

    Polos_Malha_Aberta=cplxpair(pole(sistema_malha_aberta_Gdy));

    Polos_Malha_Fechada_Alocacao_Livre=...
        cplxpair(pole(sistema_alocacao_livre_Gdy));

    Polos_Malha_Fechada_Alocacao_Parcial=...
        cplxpair(pole(sistema_alocacao_parcial_Gdy));

    Tabela_Polos=table(Polos_Malha_Aberta,...
                       Polos_Malha_Fechada_Alocacao_Parcial);

    %------------
%   15.4) Graficos para Polos e posicao 
               
    autovalores=...
        [ Polos_Malha_Aberta',...
          Polos_Malha_Fechada_Alocacao_Livre',...  
          Polos_Malha_Fechada_Alocacao_Parcial'];

    [  x_horizontal_strip,y_inferior_horizontal_strip,...
        y_superior_horizontal_strip,...
        x_sector,y_inferior_sector,y_superior_sector,...
        x_parabola,y_inferior_parabola,y_superior_parabola,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = linhas_D_regioes(...
            autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P);  

     figure
     plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
          real(Polos_Malha_Fechada_Alocacao_Livre),...
          imag(Polos_Malha_Fechada_Alocacao_Livre),'xr',...
          real(Polos_Malha_Fechada_Alocacao_Parcial),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial),'xb',...
          x_alphav,y_alphav,'-c',...
          x_betav,y_betav,'-c',...
          x_disk,y_disk,'-c',...
          x_horizontal_strip,y_inferior_horizontal_strip,'-c',...
          x_horizontal_strip,y_superior_horizontal_strip,'-c',...
          x_sector,y_inferior_sector,'-c',...
          x_sector,y_superior_sector,'-c',...
          x_parabola,y_inferior_parabola,'-c',...
          x_parabola,y_superior_parabola,'-c',...
          'LineWidth',2,'MarkerSize',8)
     xlim([-5.1 0]);
     ylim([-5.1 5.1]);

     title('Polos e localizacao Plano Complexo','FontSize',12)
     legend('Aberta','Livre','Parcial','Fronteira','FontSize',12)   
%-----------------
%  15.5) Grafico para Comparar Acao de controle a pertubacao degrau 

    [y_livre,t_livre,x_livre] = ...
        step(sistema_alocacao_livre_Gdy,0:0.01:12);
        
    u_livre(1:size(x_livre,1),1:n_u,1)=x_livre(:,:,1)*Kpf_livre';
    
    [y_parcial,t_parcial,x_parcial] = ...
        step(sistema_alocacao_parcial_Gdy,0:0.01:12);
        
    u_parcial(1:size(x_livre,1),1:n_u,1)=...
        x_parcial(:,:,1)*Kpf_parcial';

            figure
            subplot(221)
            plot(t_livre,u_livre(1:end,1),'LineWidth',2)
            title('Resposta u1(t) a d(t)','FontSize',12)
            legend('Livre','FontSize',12)
            subplot(222)
            plot(t_parcial,u_parcial(1:end,1),'LineWidth',2)
            title('Resposta u1(t) a d(t)','FontSize',12)
            legend('Parcial','FontSize',12)
            subplot(223)
            plot(t_livre,u_livre(1:end,2),'LineWidth',2)
            title('Resposta u2(t) a d(t)','FontSize',12)
            legend('Livre','FontSize',12)
            subplot(224)
            plot(t_parcial,u_parcial(1:end,2),'LineWidth',2)
            title('Resposta u2(t) a d(t)','FontSize',12)
            legend('Parcial','FontSize',12)
 %===================
 %  16) Mostrar Tabela Funcoes Parciais Malha Aberta
 
    Tabela_Funcoes_Parciais,
 %=================
 %  17) Mostrar Tabela Polos 
    Tabela_Polos,
%==================
%   18) Tabela para comparar normas dos sitemas

%   18.1) Computar normas do sistema malha aberta

    norma_H2(1,1)=norm(sistema_malha_aberta_Gdy,2);
    
    variacao_percentual_norma_H2(1,1)=0;
    
    norma_Hinf(1,1)=norm(sistema_malha_aberta_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(1,1)=0;

 %   18.2) Computar normas do sistema malha fechada alocacao livre

    norma_H2(2,1)=norm(sistema_alocacao_livre_Gdy,2);
    
    variacao_percentual_norma_H2(2,1)=...
        100*(norma_H2(2,1)-norma_H2(1,1))/norma_H2(1,1);
    
    norma_Hinf(2,1)=norm(sistema_alocacao_livre_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(2,1)=...
        100*(norma_Hinf(2,1)-norma_Hinf(1,1))/norma_Hinf(1,1);
    
%   18.3) Computar normas do sistema malha fechada alocacao parcial

    norma_H2(3,1)=norm(sistema_alocacao_parcial_Gdy,2);
    
    variacao_percentual_norma_H2(3,1)=...
        100*(norma_H2(3,1)-norma_H2(1,1))/norma_H2(1,1);
    
    norma_Hinf(3,1)=norm(sistema_alocacao_parcial_Gdz,'inf');
    
    variacao_percentual_norma_Hinf(3,1)=...
        100*(norma_Hinf(3,1)-norma_Hinf(1,1))/norma_Hinf(1,1);
    
%   18.4) Tempo de otimizacao
    Tempo_Otimizacao=[0; tempo_otimizacao_alocacao_livre;...
                       tempo_otimizacao_alocacao_parcial];
                   
    Funcoes_transferencia=["Aberta";"Alocacao Livre";"Alocacao Parcial"];
    Tabela_Normas = table(Funcoes_transferencia,...
      norma_H2,variacao_percentual_norma_H2,...
      norma_Hinf,variacao_percentual_norma_Hinf,Tempo_Otimizacao),
%---------------------------------------