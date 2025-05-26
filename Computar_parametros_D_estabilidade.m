    %A) "q_d" e "e_P"
    
    q_d=0;
    e_P=[];

    %B) "w_H" via equacoes:
    
    %; %equacao (3.58) 
    if ~isempty(Tpmin)
        w_H=pi/Tpmin;
    else
        w_H=[];
    end
    
    %C) "theta_s" via equacoes:
    
    if ~isempty(Mosmax)
        theta_s=atan( pi/log( Mosmax^( -1 ) ) ); % equacao (3.42)
    else
        theta_s=[];
    end
    
    %D) "r_d" via equacoes:
    if ~isempty(theta_s)
        if ~isempty(Trmin)
            rd1=( 1.1 + 0.125*cos( theta_s ) + 0.495*cos( theta_s )^2)/Trmin; %equacao (3.47)
        else
            rd1=[];
        end
        if ~isempty(Tdmin)
            rd2=( 1.1 - 0.416*cos( theta_s ) + 2.917*cos( theta_s )^2)/Tdmin; %equacao (3.54)
        else
            rd2=[];
        end
    else
        rd1=[];
        rd2=[];
    end
    r_d=min([rd1,rd2]);
    
    %E) "beta_v" via equacoes:
    
    if ~isempty(Trmin)
        bv1=2.2/Trmin; %equacao (3.20)
    else
        bv1=[];
    end
    
    if ~isempty(Tsmin)
        bv2=4/Tsmin; %equacao (3.21)
    else
        bv2=[];
    end
    beta_v=min([bv1,bv2]);

    %F) "alpha_v" via equacoes:
    
    %av1=2.2/Trmax; %equacaao (3.20)
    if ~isempty(Trmax)
        av1=2.2/Trmax;
    else
        av1=[];
    end
    if ~isempty(Tsmax)
    av2=4/Tsmax; %equacao (3.21)    
        if ~isempty(csmax)
            av3=( log( csc(theta_s)*csmax^(-1) ) )/Tsmax; %equacao (3.65)
        else
            av3=[];
        end
    else
       av2=[];
    end
    alpha_v=max([av1,av2,av3]);
