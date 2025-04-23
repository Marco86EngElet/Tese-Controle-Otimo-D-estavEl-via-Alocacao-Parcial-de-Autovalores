clear all,close all,clc,



%Escrevendo polos reais na diagonal principal de A_x
n_x=30;
A_x=zeros(n_x,n_x)
for i=1:n_x
    for j=1:n_x
        A_x(i,j)=cos((i+j)/5);
    end
end
 [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=...
            eig(A_x);
%Escrevendo par de polos complexo conjugado na diagonal principal de A_x
m=0;
for k=1:sqrt(n_pcc)
    for j=1:sqrt(n_pcc)
        m=m+1;
        linhas=(n_r+2*m-1):(n_r+2*m);
        colunas=linhas;
        A_x(linhas,colunas)=[-0.01*k, 0.01*j ; -0.01*j, -0.01*k]; 
    end
end

%Escrevendo matrix B_u
n_u=2;
B_u=ones(n_x,n_u);
%B_u(n_x-1:n_x,1:2)=eye(2);

%Escrevendo matriz B_d
n_d=2;
B_d=0.1*ones(n_x,n_d);

%Escrevendo matriz C_y
n_y=14;
C_y=ones(n_y,n_x);
C_y(1,n_x)=1;

%Escrevendo matriz D_y
D_y=0.01*ones(n_y,n_u);
%D_y(1,n_u)=0.001;

%Escrevendo matriz E_y
E_y=zeros(n_y,n_d);

%Escrevendo matriz C_z
n_z=15;
C_z=ones(n_z,n_x);
%C_z(1,1)=1;

%Escrevendo matriz D_z
D_z=0.01*ones(n_z,n_u);
%D_z(1,1)=0.001;

%Escrevendo matriz E_z
E_z=0.001*ones(n_z,n_d);
%E_z(1,1)=0.001;

%Verificando quantidade de estados nao controlaveis
estados_nao_controlaveis=n_x-rank(ctrb(A_x,B_u));



%Caracteristicas da região de D-estabilidade
    
alpha_v=2;
beta_v=3;
theta_s=[];
r_d=[];
q_d=[];
w_H=[];
e_P=[];

%Opcoes para otimizacao semideinida positiva no yalmip
     
Yalmip_sdpsettings =...
     sdpsettings('verbose',1,'solver','lmilab','debug',1);
c_H2=1000;
c_Hinf=1;

%Iniciar processo de otimizacao multistep
tic
Kpf_parcial=zeros(n_u,n_x);
Polos_Malha_Aberta=cplxpair(eig(A_x))

k=1;
while k<=n_x
    k,
    if imag(Polos_Malha_Aberta(k,1))==0
        
        Q=sqrt(2);
        
        [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=...
            eig(A_x+B_u*Kpf_parcial);
        
        % Defina os polos x
        polo1 = Polos_Malha_Aberta(k,1);
        
        % Calcule a diferença absoluta entre cada elemento de Lambda e
        % polos
        dif_polo1 = abs(Eigenvalues - polo1);
        
        % Encontre o índice do elemento com a menor diferença
        [~, minIndex_1] = min(dif_polo1(:));
        
        % Converta o índice linear para índices de matriz
        [row_polos, col_polos] = ind2sub(size(Eigenvalues), minIndex_1);
        
        %Encontrar as matrizes para alocacao parcial
        Lambda_p=Eigenvalues(row_polos,col_polos);
        L_p=Left_Eigenvectors(1:end,col_polos);
        
        %Matriz de controlabilidade parcial
        matriz_controlabilidade_parcial=ctrb(Lambda_p,(Q*L_p')*B_u/2);
        
        %Quantidade de estados não controlaveis sistema reduzido
        n_ola=1;
        
        estados_nao_controlaveis_sistema_reduzido=...
        n_ola-rank(matriz_controlabilidade_parcial);
        
        if estados_nao_controlaveis_sistema_reduzido==0
        
            [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz]=...
            matrizes_sistema_reduzido(Lambda_p,L_p,Q,B_u,B_d,C_y,C_z);
        
            %Executar otimizacao para alocacao parcial 
            
            controle_misto_H2_Hinf_Parcial_D_estabilidade
                    
            try
                controle_misto_H2_Hinf_Parcial_D_estabilidade
                if ~isnan(KD_parcial)
                    if ~isinf(KD_parcial)
                        tilde_Lp=Q*L_p'/2;
                        Kpf_parcial=Kpf_parcial+KD_parcial*tilde_Lp;
                    end
                end
            catch
                
            end
        end
        k=k+1;
    else
        k
        Q = [1 1;1i -1i];
        
        [Right_Eigenvectors,Eigenvalues,Left_Eigenvectors]=...
            eig(A_x+B_u*Kpf_parcial);
        Polos=diag(Eigenvalues);
        
        % Defina os polos x
        polo1 = Polos_Malha_Aberta(k,1);
        k=k+1,
        polo2 = Polos_Malha_Aberta(k,1);
        
        % Calcule a diferença absoluta entre cada elemento de Lambda e
        % polos
        dif_polo1 = abs(Polos - polo1);
        dif_polo2 = abs(Polos - polo2);
        
        % Encontre o índice do elemento com a menor diferença
        [~, minIndex_1] = min(dif_polo1(:));
        [~, minIndex_2] = min(dif_polo2(:));
        
        %Agrupar os indices 
        min_row=min(minIndex_1,minIndex_2);
        max_row=max(minIndex_1,minIndex_2)
        row_polos=min_row:max_row;
        col_polos=row_polos;
        
        %Encontrar as matrizes para alocacao parcial
        Lambda_p=Eigenvalues(row_polos,col_polos);
        L_p=Left_Eigenvectors(1:end,col_polos);
        
        n_ola=2;
        
        matriz_controlabilidade_parcial=ctrb(Lambda_p,(Q*L_p')*B_u/2);
        estados_nao_controlaveis_sistema_reduzido=...
        n_ola-rank(matriz_controlabilidade_parcial);
        
        if estados_nao_controlaveis_sistema_reduzido==0
            
            %Quantidade de estados não controlaveis sistema reduzido
            [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz]=...
            matrizes_sistema_reduzido(Lambda_p,L_p,Q,B_u,B_d,C_y,C_z);

            %Executar otimizacao para alocacao parcial 
            try
                controle_misto_H2_Hinf_Parcial_D_estabilidade
                if ~isnan(KD_parcial)
                    if ~isinf(KD_parcial)
                        tilde_Lp=Q*L_p'/2;
                        Kpf_parcial=Kpf_parcial+KD_parcial*tilde_Lp;
                    end
                end
            catch
                
            end
        end
        k=k+1;
    end
end
tempo_otimizacao=toc
%==========================================================================    
%   14) Sistema malha fechada via controle otimo alocacao parcial
    
    sistema_malha_aberta_y=ss(A_x,B_d,C_y,E_y);
    
    sistema_malha_aberta_z=ss(A_x,B_d,C_z,E_z);

    sistema_alocacao_parcial_y=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_y+D_y*Kpf_parcial,E_y);
    
    sistema_alocacao_parcial_z=...
        ss(A_x+B_u*Kpf_parcial,B_d,C_z+D_z*Kpf_parcial,E_z);
    
%==========================================================================
%   17) Computar normas do sistema malha aberta

    norma_H2(1,1)=norm(sistema_malha_aberta_y,2);
    
    variacao_percentual_norma_H2(1,1)=0;
    
    norma_Hinf(1,1)=norm(sistema_malha_aberta_z,'inf');
    
    variacao_percentual_norma_Hinf(1,1)=0;

%   19) Computar normas do sistema malha fechada alocacao parcial

    norma_H2(2,1)=norm(sistema_alocacao_parcial_y,2);
    
    variacao_percentual_norma_H2(2,1)=...
        100*(norma_H2(2,1)-norma_H2(1,1))/norma_H2(1,1);
    
    norma_Hinf(2,1)=norm(sistema_alocacao_parcial_z,'inf');
    
    variacao_percentual_norma_Hinf(2,1)=...
        100*(norma_Hinf(2,1)-norma_Hinf(1,1))/norma_Hinf(1,1);
    
%   20) Tabela para comparar normas

    Funcoes_transferencia=["Aberta";"Alocacao Parcial MultiStep"];
    Tabela_Normas = table(Funcoes_transferencia,...
        norma_H2,variacao_percentual_norma_H2,...
        norma_Hinf,variacao_percentual_norma_Hinf),
%==========================================================================
%   26) Grafico para Polos e posicao 

     Polos_Malha_Fechada_Alocacao_Parcial_MultiStep=...
        cplxpair(pole(sistema_alocacao_parcial_y));
    
    autovalores=[ Polos_Malha_Aberta',...
                  Polos_Malha_Fechada_Alocacao_Parcial_MultiStep'];

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
          real(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),...
          imag(Polos_Malha_Fechada_Alocacao_Parcial_MultiStep),'xb',...
          x_alphav,y_alphav,'-k',...
          x_betav,y_betav,'-k',...
          x_sector,y_inferior_sector,'-k',...
          x_sector,y_superior_sector,'-k')    