%29.5) Obter a Matriz de Controlabilidade Parcial (Matriz de 
    %   Controlabilidade do SREAPE) 

    controlability_matrix = ctrb(Lambda_p, (Q_p * L_p') * B_u / 2);

    %29.6) Calcular numeros de estados nao controlaveis do SREAPE

    number_of_non_controlable_states = 2-rank(controlability_matrix);

    %29.7) Tentar executar o Algoritmo 7
    if number_of_non_controlable_states == 0 

        %29.8) Computar matrizes onstante oriundas das 
        %   transformacoes de similaridades

        Tilde_Lambda_p = Q_p*Lambda_p*Q_p'/2;
        Tilde_Bu = Q_p*L_p'*B_u/2;
        Tilde_Bd = Q_p*L_p'*B_d/2;

        if ~isempty(C_y)

            Tilde_Cy = C_y*L_p*inv(L_p'*L_p)*Q_p'/2;

        else
            Tilde_Cy=[];
        end
        if ~isempty(C_z)

            Tilde_Cz = C_z*L_p*inv(L_p'*L_p)*Q_p'/2;

        else

            Tilde_Cz=[];

        end

        %29.9) Construir variaveis de decisao

        Tilde_X  = sdpvar(n_ola,n_ola,'symmetric');
        Tilde_W  = sdpvar(n_u,n_ola,'full');

        %29.10) Iniciar LMIs para Algoritmo 7 para alocar o 
        %   j-enesimo par de polos complexos

        set_LMIs_parcial = Tilde_X>=eps*eye(n_ola);

        %29.11) Construir LMI para "Re(s)<=-alpha_v" (Regiao Faixa
        %   Vertical)
        if ~isempty(alpha_v)

            Tilde_Ax_Valpha = Tilde_Lambda_p+alpha_v*eye(n_ola);

            set_LMIs_parcial = [ set_LMIs_parcial,...
                Tilde_Ax_Valpha*Tilde_X+Tilde_X*Tilde_Ax_Valpha'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola)  ];
        end

        %29.12) Construir LMI para "Re(s)>=-beta_v" (Regiao Faixa
        %   Vertical)
        if ~isempty(beta_v)

            Tilde_Ax_Vbeta = -Tilde_Lambda_p-beta_v*eye(n_ola);

            set_LMIs_parcial = [ set_LMIs_parcial,...
                Tilde_Ax_Vbeta*Tilde_X+Tilde_X*Tilde_Ax_Vbeta'-...
                Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu'<=-eps*eye(n_ola) ];
        end

        %29.13) Construir LMI para "abs(s)<=-r_d" (Regiao Disco)
        if ~isempty(r_d)

            Tilde_Ax_Dqr = Tilde_Lambda_p+q_d*eye(n_ola); 

            set_LMIs_parcial = [ set_LMIs_parcial,...
                [ -r_d*Tilde_X,...
                  Tilde_Ax_Dqr*Tilde_X+Tilde_Bu*Tilde_W;...
                  Tilde_X*Tilde_Ax_Dqr'+Tilde_W'*Tilde_Bu',...
                  -r_d*Tilde_X ]...
                  <=-eps*eye(2*n_ola)    ];
        end

        %29.14) Construir LMI para "-theta_s<=angle(s)<=theta_s"
        % (Regiao Setor)
        if ~isempty(theta_s)

            Tilde_Axsin = Tilde_Lambda_p*sin(theta_s);

            Tilde_Axcos = Tilde_Lambda_p*cos(theta_s);

            Tilde_Busin = Tilde_Bu*sin(theta_s);

            Tilde_Bucos = Tilde_Bu*cos(theta_s);

            set_LMIs_parcial =[ set_LMIs_parcial,...
                [ Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin',...
                  Tilde_Axcos*Tilde_X-Tilde_X*Tilde_Axcos'+...
                  Tilde_Bucos*Tilde_W-Tilde_W'*Tilde_Bucos';...
                  Tilde_X*Tilde_Axcos'-Tilde_Axcos*Tilde_X+...
                  Tilde_W'*Tilde_Bucos'-Tilde_Bucos*Tilde_W,...
                  Tilde_Axsin*Tilde_X+Tilde_X*Tilde_Axsin'+...
                  Tilde_Busin*Tilde_W+Tilde_W'*Tilde_Busin'...
                ]<=-eps*eye(2*n_ola) ];
        end

        %29.15) Construir LMI para "-w_H<=imag(s)<=w_H" (Regiao 
        %   Faixa Horizontal)
        if ~isempty(w_H)
            set_LMIs_parcial =[ set_LMIs_parcial,...
                [ -w_H*Tilde_X,...
                  Tilde_X*Tilde_Lambda_p'-Tilde_Lambda_p*Tilde_X+...
                  Tilde_W'*Tilde_Bu'-Tilde_Bu*Tilde_W;...
                  Tilde_Lambda_p*Tilde_X-Tilde_X*Tilde_Lambda_p'+...
                  Tilde_Bu*Tilde_W-Tilde_W'*Tilde_Bu',...
                  -w_H*Tilde_X...
                ]<=-eps*eye(2*n_ola) ];
        end 

        %29.16) Construir LMI para Control otimo H_2 
        if ~isempty(c_H2)

            Z = sdpvar(n_y,n_y,'symmetric');

            rho = sdpvar(1,1,'symmetric');   

            set_LMIs_parcial = [ set_LMIs_parcial, ...
                Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu'+...
                Tilde_Bd*Tilde_Bd'<=-eps*eye(n_ola) ];

            set_LMIs_parcial = [ set_LMIs_parcial,...
                        Z>=eps*eye(n_y),...
                        trace(Z)<=rho,...
                        rho>=eps ];

            set_LMIs_parcial = [set_LMIs_parcial,...
                [ -Z, Tilde_Cy*Tilde_X+D_y*Tilde_W;...
                    Tilde_X'*Tilde_Cy'+Tilde_W'*D_y', -Tilde_X ]<=...
                    -eps*eye(n_y+n_ola,n_y+n_ola) ]        
        end

        %29.17) Construir LMI para Controle otimo H_infinito
        if ~isempty(c_Hinf)

            gamma = sdpvar(1,1,'symmetric');   

            set_LMIs_parcial = [ set_LMIs_parcial,...
                gamma>=eps,...
                [   Tilde_Lambda_p*Tilde_X+Tilde_X*Tilde_Lambda_p'+...
                    Tilde_Bu*Tilde_W+Tilde_W'*Tilde_Bu',...
                    Tilde_Bd, Tilde_X*Tilde_Cz'+Tilde_W'*D_z';...
                    Tilde_Bd', -gamma*eye(n_d), E_z';...
                    Tilde_Cz*Tilde_X+D_z*Tilde_W,...
                    E_z, -gamma*eye(n_z)...
                ]<=-eps*eye(n_ola+n_d+n_z)...
            ];
        end

        else %Caso SREAPE seja incontrolavel
            disp('SREAPE incontrolavel')
        end        