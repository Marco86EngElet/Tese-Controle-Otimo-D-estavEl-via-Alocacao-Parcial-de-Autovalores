%34.1) Sinal (y(t))^2 gerado.
    
    [y_Gcdys_impulse_parcial,t_Gcdys_impulse_parcial,x_Gcdys_impulse_parcial]=...
        impulse(Gcdys_parcial,tempo_simulacao_impulso);
    
    y2_Gcdys_impulse_parcial=y_Gcdys_impulse_parcial.*y_Gcdys_impulse_parcial;
    
    %34.2) Sinal u(t) gerado.
    
    u_Gcdys_impulse_parcial(1:size(x_Gcdys_impulse_parcial,1),1:n_u,1)=...
        x_Gcdys_impulse_parcial(:,:,1)*Kpf_parcial';

    u2_Gcdys_impulse_parcial=u_Gcdys_impulse_parcial.*u_Gcdys_impulse_parcial;