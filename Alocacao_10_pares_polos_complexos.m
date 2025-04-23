clear all,close all,clc,

%   1) Coeficiente de armotecimento, Frequencia Natural e Ganho das funcoes
%       funcoes parciais
%       1.1) Coeficiente de amortecimento
zeta_c(1,1)=0.1; %Funcao 1
zeta_c(2,1)=0.12; %Funcao 2
zeta_c(3,1)=0.15; %Funcao 3
zeta_c(4,1)=0.2; %Funcao 4
zeta_c(5,1)=0.3; %Funcao 5

%       1.2) Frequencia nao natural
w_c(1,1)=0.1; %Funcao 1
w_c(2,1)=0.12; %Funcao 2
w_c(3,1)=0.15; %Funcao 3
w_c(4,1)=0.2; %Funcao 4
w_c(5,1)=0.3;%Funcao 5

%   2)  Garacteristicas dos polos de cada fracao parcial
%
%       2.1) Parte real associada wn(i,1) e zeta(i,1)
for c=1:5
    parte_real(c,1)=w_c(c)*zeta_c(c,1);
end

%       2.2) Parte imaginaria associada wn(i,1) e zeta(i,1)
for c=1:5
    parte_imaginaria(c,1)=w_c(c,1)*(sqrt(1-zeta_c(c,1)^2));
end

%       2.3) Angulo associado em graus a zeta(i,1)
for c=1:5
    angulo_degraus(c,1)=rad2deg(acos(zeta_c(c,1)));
end
%   3)  Medidas transitorias associadas aos polos
%
%       3.1) Criterio de amortecimento
    csmax=0.02;
%       3.2) Medidas transitorias
    %mos -> maximo sobresinal ('overshoot')
    %ts -> tempo de acomodacao ('settling time')
    %td -> tempo de atraso ('delay time')
    %tr -> tempo de subida ('rise time')
    %tp -> tempo de pico ('peak time')
for c=1:5    
    [mos(c,1),ts(c,1),td(c,1),tr(c,1),tp(c,1)] =...
        medidas_SISO_2nd(zeta_c(c,1),w_c(c,1),csmax);
end

%   4)  Criar funcoes parciais
for c=1:5
G_Ycij(c)=...
Criar_Funcao_Transferencia_SISO_2nd(1,zeta_c(c,1),w_c(c,1)); 
end

%   5)  Tabela para descrever caracteristicas transitoria e dos polos das
%           funcoes parciais

Tabela_Funcoes_Parciais=...
table(zeta_c,w_c,parte_real,parte_imaginaria,angulo_degraus,...
  mos,ts,td,tr,tp);

%   6) Computando as funcoes de saidas 

%       6.1) Entrada 'd1(t)' e saida 'y1(t)'
G_Y11=G_Ycij(1)+G_Ycij(2)+G_Ycij(3)+G_Ycij(4)+G_Ycij(5);

%       6.4) Entrada 'd1(t)' e saida 'z1(t)'
G_Z11=G_Ycij(1)+G_Ycij(2)-G_Ycij(3)+G_Ycij(4)+G_Ycij(5)+0.01;

%   7)  Computando modelo em espaco de estados

%       7.1) Obter espaco de estados para 'G_d1y' 
Gd1y1_tf_ss = ss(G_Y11);

%       7.2) Obter espaco de estados para 'G_d1z'
Gd1z1_tf_ss = ss(G_Z11);
A_x=Gd1y1_tf_ss.A;

%       7.4) Matriz B_d para modelo de estado de estados, B_d=[B1,B2,B3]
B_d=Gd1y1_tf_ss.B;

%       7.5) Matriz C_y para modelo de estado de estados, C_y=[C1;C2;C3]
C_y=Gd1y1_tf_ss.C; 

%       7.6) Matriz E_y para modelo de estado de estados, E_y=[E1,E2,E3]
E_y=blkdiag(Gd1y1_tf_ss.D);

%       7.7) Matriz C_z para modelo de estado de estados, C_z=[C1;C2;C3]
C_z=Gd1z1_tf_ss.C;

%       7.8) Matriz E_z para modelo de estado de estados, E_y=[E1,E2,E3]
E_z=blkdiag(Gd1z1_tf_ss.D); 

%       7.9) Objeto para sistema em malha aberta para Gdy        
sistema_malha_aberta_Gdy=ss(A_x,B_d,C_y,E_y);

%       7.10) Objeto para sistema em malha aberta para Gdz  
sistema_malha_aberta_Gdz=ss(A_x,B_d,C_z,E_z);

%   8) Obter dimensoes dos vetores de estados, entradas e saidas
n_x=size(A_x,1);
n_y=size(C_y,1);
n_z=size(C_z,1);
n_d=size(B_d,2);    
n_u=2;
n_ola=n_u;
%   9) Matrizes associadas as acoes de controladores no modelo em espaco de 
%       estados
%   9.1) Matriz B_u
B_u=zeros(n_x,n_u);
B_u(2,1)=1;
B_u(3,1)=1;
B_u(2,2)=-1;
B_u(3,2)=-1;

%       9.2)  Matriz D_y
D_y=0.001*ones(n_y,n_u);

%       9.3)  Matriz D_z
D_z=D_y;

%   10) Caracteristicas da regiao de D-estabilidade

%       10.1) Parte Imaginaria entre '-w_H' e "+w_H"
Tpmin=2.5;
w_H=pi/Tpmin;

%       10.2) Angulo dos polos entre '-theta_s' e '+theta_s'
mosmax=0.2;
theta_s=atan(pi/log(mosmax^(-1)));

%       10.3) Raio disco
Trmin=1.8;
r_d=(1.1+0.125*cos(theta_s)+0.495*cos(theta_s)^2)/Trmin;

%       10.4) Parte real inferior a '-alpha_v'
Trmax=6;
Tsmax=8.8;
av1=2.2/Trmax;
av2=4/Tsmax;
av3=(log(csc(theta_s)*csmax^(-1)))/Tsmax;
alpha_v=max([av1,av2,av3]);

%       10.5) Parte real superior a '-beta_v'
Tsmin=5;
bv1=2.2/Trmin;
bv2=4/Tsmin;
beta_v=min([bv1,bv2]);

%       10.6) Centro Disco
q_d=0;

%       10.7) Fator de amortecimento Parabola de Estabilidade
e_P=[];

%   11) Configuracoes no yalmip
%11.1) Configuracoes para otimizacao semidefinida no yalmip
Yalmip_sdpsettings =...
     sdpsettings('verbose',1,'solver','lmilab','debug',1);

%11.2) Pesos da funcao custo de otimizacao 
c_H2=1;
c_Hinf=2;

%   12) Computar matriz de retroalimentacao com alocacao livre
%       12.1) Matriz de Controlabilidade Total
matriz_controlabilidade_total=ctrb(A_x,B_u);

%       12.2) Quantidade de estados nao controlaveis para sistema original
estados_nao_controlaveis_sistema_original=...
n_x-rank(matriz_controlabilidade_total);

%       12.3) Executar controle misto H2/Hinf classico
if estados_nao_controlaveis_sistema_original==0
    tic
    controle_misto_H2_Hinf_D_estabilidade    
    tempo_otimizacao_alocacao_livre=toc;
else
   disp('sistema nao controlavel') 
end

%   13) Sistema malha fechada via controle otimo alocacao livre

%   13.1) Sistema d(t) para y(t) 
sistema_alocacao_livre_Gdy=...
ss(A_x+B_u*Kpf_livre,B_d,C_y+D_y*Kpf_livre,E_y);

%   13.2) Sistema d(t) para z(t) 
sistema_alocacao_livre_Gdz=...
ss(A_x+B_u*Kpf_livre,B_d,C_z+D_z*Kpf_livre,E_z); 

%   14) Computar matriz de retroalimentacao com alocacao parcial
%       14.1) Matriz de transformacao
Q = [1 1;1i -1i];

%       14.2) Iniciando matriz de retroalimentacao
Kpf_parcial=zeros(n_u,n_x);

%       14.3) Obter polos do sistema em malha
Polos_Malha_Aberta=cplxpair(pole(sistema_malha_aberta_Gdy));

%       14.4) Iniciar contagem de tempo de otimizacao para alocacao parcial
tic

%       14.5) Estrutura de contagem para alocacao parcial multi-etapa
for j=1:5
%           14.6) Obter Matriz dos Autovalores e Autovetores
    [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=...
        eig(A_x+B_u*Kpf_parcial);

%           14.7) Defina os polos x
    polo1 = Polos_Malha_Aberta(2*j-1,1);
    polo2 = Polos_Malha_Aberta(2*j,1);

%           14.8) Calcule a diferenca absoluta entre cada elemento de 
%           Lambda e polos
    dif_polo1 = abs(Eigenvalues - polo1);
    dif_polo2 = abs(Eigenvalues - polo2);

%           14.9) Encontre o indice do elemento com a menor diferenca
    [~, minIndex_1] = min(dif_polo1(:));
    [~, minIndex_2] = min(dif_polo2(:));

%           14.10) Converta o indice linear para indices de matriz
    [row_1, col_1] = ind2sub(size(Eigenvalues), minIndex_1);
    [row_2, col_2] = ind2sub(size(Eigenvalues), minIndex_2);

    %           14.11) Agrupar os indices 
    row_polos=min(row_1,row_2):max(row_1,row_2);
    col_polos=min(col_1,col_2):max(col_1,col_2);

    %           14.12) Encontrar as matrizes para alocacao parcial
    Lambda_p=Eigenvalues(row_polos,col_polos);
    L_p=Left_Eigenvectors(1:end,col_polos);

    %           14.13) Obter matriz de controlabidade parcial
    matriz_controlabilidade_parcial=ctrb(Lambda_p,(Q*L_p')*B_u/2);

    %           14.14) Calcular quantidade de estados nao controlaveis sistema 
%           reduzido
    estados_nao_controlaveis_sistema_reduzido=...
    n_ola-rank(matriz_controlabilidade_parcial);

%           14.15) Executar otimizacao para alocacao parcial
    if estados_nao_controlaveis_sistema_reduzido==0

        %           14.16) Obter matrizes do sistema reduzido
        [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz]=...
        matrizes_sistema_reduzido(Lambda_p,L_p,Q,B_u,B_d,C_y,C_z);

    % 14.17) Executar otimizacao para alocacao parcial 
        controle_misto_H2_Hinf_Parcial_D_estabilidade

        % 14.18) Obter matriz de retroalimentacao
        Kpf_parcial=Kpf_parcial+KD_parcial*(Q*L_p')/2;
    else
        disp('polo nao controlavel detectado')
    end
    
end

%       14.19) Encerrar tempo de otimizacao multietapa        
tempo_otimizacao_alocacao_parcial_multistep=toc;

%   15) Sistema malha fechada via controle otimo alocacao parcial

%   15.1) Sistema d(t) para y(t)
sistema_alocacao_parcial_y=...
ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);

%   15.2) Sistema d(t) para z(t)
sistema_alocacao_parcial_z=...
ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z);   

%   16) Mostrar Tabela Funcoes Parciais Malha Aberta
Tabela_Funcoes_Parciais,

%   17) Mostrar Tabela Polos 
%       17.1) Polos malha fechada via alocacao livre
    Polos_Malha_Fechada_Alocacao_Livre=...
        cplxpair(pole(sistema_alocacao_livre_Gdy));
    
%       17.2) Polos malha fechada alocacao parcial multistep   
    Polos_Malha_Fechada_Alocacao_Parcial_MultiStep=...
        cplxpair(pole(sistema_alocacao_parcial_y));
        Tabela_Polos=table(Polos_Malha_Aberta,...
               Polos_Malha_Fechada_Alocacao_Livre,...
               Polos_Malha_Fechada_Alocacao_Parcial_MultiStep);
           
%       17.3) Mostrar tabela                   
    Tabela_Polos,
    
%   17) Computar normas do sistema malha aberta
%       17.1) norma H2 para Gdy em malha aberta
    norma_H2(1,1)=norm(sistema_malha_aberta_Gdy,2);
    
%       17.2) variacao norma H2 para Gdy em malha aberta            
    variacao_percentual_norma_H2(1,1)=0;
    
%       17.3) norma Hinf para Gdy em malha aberta            
    norma_Hinf(1,1)=norm(sistema_malha_aberta_Gdz,'inf');
    
%       17.4) variacao norma Hinf para Gdinf em malha aberta            
    variacao_percentual_norma_Hinf(1,1)=0;
    
%       17.5) Funcao_Custo em malha aberta            
    Funcao_custo(1,1)=c_H2*norma_H2(1,1)+c_Hinf*norma_Hinf(1,1);
    
%       17.6) variacao norma Hinf para Gdinf em malha aberta            
    variacao_percentual_Funcao(1,1)=0;
    
%   18) Computar normas do sistema malha fechada alocacao livre
%       18.1) norma H2 para Gdy em malha aberta
norma_H2(2,1)=norm(sistema_alocacao_livre_Gdy,2);

%       18.2) variacao norma H2 para Gdy em malha aberta            
variacao_percentual_norma_H2(2,1)=...
    100*(norma_H2(2,1)-norma_H2(1,1))/norma_H2(1,1);

%       18.3) norma Hinf para Gdy em malha aberta            
norma_Hinf(2,1)=norm(sistema_alocacao_livre_Gdz,'inf');

%       18.4) variacao norma Hinf para Gdinf em malha aberta            
variacao_percentual_norma_Hinf(2,1)=...
    100*(norma_Hinf(2,1)-norma_Hinf(1,1))/norma_Hinf(1,1);

%       18.5) Funcao_Custo alocacao livre            
    Funcao_custo(2,1)=c_H2*norma_H2(2,1)+c_Hinf*norma_Hinf(2,1);
    
%       18.6) variacao norma Hinf para alocacao livre            
    variacao_percentual_Funcao(2,1)=0;
    
%   19) Computar normas do sistema malha fechada alocacao parcial
%       19.1) norma H2 para Gdy em malha aberta
norma_H2(3,1)=norm(sistema_alocacao_parcial_y,2);

%       19.2) variacao norma H2 para Gdy em malha aberta            
variacao_percentual_norma_H2(3,1)=...
    100*(norma_H2(3,1)-norma_H2(1,1))/norma_H2(1,1);

%       19.3) norma Hinf para Gdy em malha aberta            
norma_Hinf(3,1)=norm(sistema_alocacao_parcial_z,'inf');

%       19.4) variacao norma Hinf para Gdinf em malha aberta            
variacao_percentual_norma_Hinf(3,1)=...
    100*(norma_Hinf(3,1)-norma_Hinf(1,1))/norma_Hinf(1,1);

%       19.5) Funcao_Custo alocacao livre            
    Funcao_custo(3,1)=c_H2*norma_H2(3,1)+c_Hinf*norma_Hinf(3,1);
    
%       19.6) variacao norma Hinf para alocacao livre            
    variacao_percentual_Funcao(3,1)=0;
    
%   20) Computar Tempo de Otimizacao
%       20.1) Tempo de Oitmizacao malha aberta
    Tempo_Otimizacao(1,1)=0;
    
%       20.2) Tempo de Oitmizacao multistep
    Tempo_Otimizacao(2,1)=...
        tempo_otimizacao_alocacao_livre; 
    
%       20.3) Tempo de Oitmizacao multistep
    Tempo_Otimizacao(3,1)=...
        tempo_otimizacao_alocacao_parcial_multistep; 

    
%   21) Tabela para comparar normas
Funcoes_transferencia=["Aberta";"Alocacao Livre";...
                    "Alocacao Parcial MultiStep"];
Tabela_Normas = table(Funcoes_transferencia,...
norma_H2,variacao_percentual_norma_H2,...
norma_Hinf,variacao_percentual_norma_Hinf,...
Tempo_Otimizacao),

%   22) Tabela Tempo de otimizacao
diferenca_percentual_tempo_otimizacao=...
100*(tempo_otimizacao_alocacao_parcial_multistep-...
     tempo_otimizacao_alocacao_livre)/...    
     tempo_otimizacao_alocacao_livre;
Tabela_Tempo_Otimizacao=table(...
tempo_otimizacao_alocacao_livre,...
tempo_otimizacao_alocacao_parcial_multistep,...
diferenca_percentual_tempo_otimizacao), 

%   23) Tabela Integral Numerico do Valor Absoluto da Acao de Controle
%      23.1) Resposta Transitoria alocacao livre
[y_livre,t_livre,x_livre] = ...
    step(sistema_alocacao_livre_Gdy,0:0.01:14);

%      23.2) acao de controle alocacao livre         
u_livre(1:size(x_livre,1),1:n_u,1)=x_livre(:,:,1)*Kpf_livre';

%      23.3) Integrais acoes de controle alocacao livre
Integral_Alocacao_Livre(1,1) = trapz(t_livre,...
                                 abs(u_livre(1:end,1,1)));
Integral_Alocacao_Livre(2,1) = trapz(t_livre,...
                                 abs(u_livre(1:end,2,1)));
Integral_Alocacao_Livre(3,1) = sum(Integral_Alocacao_Livre); 

%      23.4) Resposta Transitoria alocacao parcial multstep                                 
[y_parcial,t_parcial,x_parcial] = ...
    step(sistema_alocacao_parcial_y,0:0.01:14);

%      23.5) acao de controle alocacao parcial multstep 
u_parcial(1:size(x_livre,1),1:n_u,1)=...
    x_parcial(:,:,1)*Kpf_parcial';

%      23.6) Integrais acoes de controle alocacao parcial multistep
Integral_Alocacao_Parcial(1,1) = trapz(t_parcial,...
                                 abs(u_parcial(1:end,1,1)));
Integral_Alocacao_Parcial(2,1) = trapz(t_parcial,...
                                 abs(u_parcial(1:end,2,1)));
Integral_Alocacao_Parcial(3,1) = sum(Integral_Alocacao_Parcial);

%       23.7)  Discrepancia acoes de controle
Discrepencia_Integral_Porcento=...
    Integral_Alocacao_Parcial-Integral_Alocacao_Livre;
Discrepencia_Integral_Porcento=...
    Discrepencia_Integral_Porcento./Integral_Alocacao_Livre;
Discrepencia_Integral_Porcento=100*Discrepencia_Integral_Porcento;

%       23.8) Tabela
Variavel_saida =["u1(t)";"u2(t)";"u1(t)+u2(t)"];
Tabela_Integral_vetor_u=...
    table(Variavel_saida,...
        Integral_Alocacao_Livre,...
        Integral_Alocacao_Parcial,...
        Discrepencia_Integral_Porcento),

%   24) Grafico Resposta Degrau

%24.1) Configuracoes para graficos
time_opt=timeoptions;
time_opt.InputLabels.FontSize=12;
time_opt.OutputLabels.FontSize=12;
time_opt.XLabel.FontSize=12;
time_opt.YLabel.FontSize=12;
time_opt.TickLabel.FontSize=12;

%       24.2) Grafico resposta y(t)    
figure
subplot(121)
stepplot(sistema_alocacao_livre_Gdy,'-r',...
     sistema_alocacao_parcial_y,'-.b',...
     sistema_malha_aberta_Gdy,'-k',time_opt)
title('Resposta degrau y(t)','FontSize',12)
legend('Livre','Parcial','Aberta','FontSize',10)
subplot(122)
stepplot(sistema_alocacao_livre_Gdz,'-r',...
     sistema_alocacao_parcial_z,'-.b',...
     sistema_malha_aberta_Gdz,'-k',time_opt)
title('Resposta degrau z(t)','FontSize',12)
legend('Livre','Parcial','Aberta','FontSize',10)

%       24.3) Grafico resposta y(t) para sistemas malha fechada    
figure
subplot(121)
stepplot(sistema_alocacao_livre_Gdy,'-r',...
     sistema_alocacao_parcial_y,'-.b',time_opt)
title('Resposta degrau y(t)','FontSize',12)
legend('Livre','Parcial','FontSize',10)

%       24.4) Grafico resposta z(t) para sistemas malha fechada        
subplot(122)
stepplot(sistema_alocacao_livre_Gdz,'-r',...
     sistema_alocacao_parcial_z,'-.b',time_opt)
title('Resposta degrau z(t)','FontSize',12)
legend('Livre','Parcial','FontSize',10)

%   25) Grafico Bode
%       25.0) Configuraes para Bode
bode_options=bodeoptions;
bode_options.PhaseVisible='off';
bode_options.InputLabels.FontSize=12;
bode_options.OutputLabels.FontSize=12;
bode_options.Title.FontSize=12;
bode_options.XLabel.FontSize=12;
bode_options.YLabel.FontSize=12;
bode_options.TickLabel.FontSize=12;

%       25.1) Grafico resposta Gdy(jw) 
figure
subplot(121)
bode(sistema_alocacao_livre_Gdy,'-r',...
     sistema_alocacao_parcial_y,'-b',...
     sistema_malha_aberta_Gdy,bode_options,'-k',bode_options)
title('Bodemag Gdy','FontSize',12)
legend('Livre','Parcial','Aberta','FontSize',10)
subplot(122)
bodemag(sistema_alocacao_livre_Gdy,'-r',...
     sistema_alocacao_parcial_y,'-b',bode_options)
title('Bodemag Gdy','FontSize',12)
legend('Livre','Parcial','FontSize',10)

%       25.2) Grafico resposta z(t)    
figure
subplot(121)
bodemag(sistema_alocacao_livre_Gdz,'-r',...
     sistema_alocacao_parcial_z,'-b',...
     sistema_malha_aberta_Gdz,'-k',bode_options)
title('Bodemag Gdz','FontSize',12)
legend('Livre','Parcial','Aberta','FontSize',10)
subplot(122)
bodemag(sistema_alocacao_livre_Gdz,'-r',...
     sistema_alocacao_parcial_z,'-b',bode_options)
title('Bodemag Gdz','FontSize',12)
legend('Livre','Parcial','FontSize',10)

%   26) Figura acoes de controle
figure
%       26.1) Resposta u1(t) para entrada d1(t) degrau    
subplot(211)
plot(t_livre,u_livre(1:end,1,1),'-r',...
     t_parcial,u_parcial(1:end,1,1),'-b')
title('Resposta u1(t) ao degrau d1(t)','FontSize',12)
legend( 'livre','parcial','FontSize',10)

%       26.2) Resposta u2(t) para entrada d1(t) degrau            
subplot(212)
plot(t_livre,u_livre(1:end,2,1),'-r',...
     t_parcial,u_parcial(1:end,2,1),'-b')
title('Resposta u2(t) ao degrau d1(t)','FontSize',12)
legend( 'livre','parcial','FontSize',10)

%   27) Grafico para Polos e posicao 

%       27.1) Autovalores (Polos) a serem inseridos no grafico
autovalores=[ Polos_Malha_Aberta',...
              Polos_Malha_Fechada_Alocacao_Livre',...  
              Polos_Malha_Fechada_Alocacao_Parcial_MultiStep'];

%       27.2) Computando curvas de fronteira das regioes    
[  x_horizontal_strip,y_inferior_horizontal_strip,...
        y_superior_horizontal_strip,...
        x_sector,y_inferior_sector,y_superior_sector,...
        x_parabola,y_inferior_parabola,y_superior_parabola,...
        x_disk,y_disk,...
        x_alphav,y_alphav,...
        x_betav,y_betav ] = linhas_D_regioes(...
        autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P); 

%       27.3) Grafico            
    figure
    subplot(121)
    plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
         real(Polos_Malha_Fechada_Alocacao_Livre),...
         imag(Polos_Malha_Fechada_Alocacao_Livre),'xr',...
         real(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),...
         imag(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),...
         'xb',x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-c',...
         x_horizontal_strip,y_superior_horizontal_strip,'-c',...
         'LineWidth',2,'MarkerSize',8);
xlim([-0.81 0]);
ylim([-1.3 1.3]);
title('Polos e localizacao Plano Complexo','FontSize',12)
legend('Aberta','Livre','Parcial','Fronteira','FontSize',10)
    subplot(122)
    plot(real(Polos_Malha_Aberta),imag(Polos_Malha_Aberta),'ok',...
         real(Polos_Malha_Fechada_Alocacao_Livre),...
         imag(Polos_Malha_Fechada_Alocacao_Livre),'xr',...
         real(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),...
         imag(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),...
         'xb',x_alphav,y_alphav,'-c',...
         x_betav,y_betav,'-c',...
         x_sector,y_inferior_sector,'-c',...
         x_sector,y_superior_sector,'-c',...
         x_disk,y_disk,'-c',...
         x_horizontal_strip,y_inferior_horizontal_strip,'-c',...
         x_horizontal_strip,y_superior_horizontal_strip,'-c',...
         'LineWidth',2,'MarkerSize',8);
xlim([-0.46 -0.455]);
ylim([-0.5295 -0.528]);
title('Polos e localizacao Plano Complexo','FontSize',12)
legend('Aberta','Livre','Parcial','Fronteira','FontSize',10)


