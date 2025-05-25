%14.1) Matriz de Controlabilidade Total
    
    matriz_controlabilidade_total=ctrb(A_x,B_u);

    %14.2) Quantidade de estados nao controlaveis para sistema original
    
    estados_nao_controlaveis_sistema_original=...
        n_x-rank(matriz_controlabilidade_total);