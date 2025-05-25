%11.1) Sinal y(t), x(t) e o 't' associado a d(t) impulso
    
    [y_Canon_Godys_impulse,t_Canon_Godys_impulse,x_Canon_Godys_impulse]=...
        impulse(Canon_Godys,tempo_simulacao_impulso);
    
    %11.2) Sinal (y(t))^2 associado a d(t) impulso
    
    y2_Canon_Godys_impulse=y_Canon_Godys_impulse.*y_Canon_Godys_impulse;