function [Tilde_Ax,Tilde_Bu,Tilde_Bd,Tilde_Cy,Tilde_Cz] =...
    matrizes_sistema_reduzido(Lambda,L_j,Q,B_u,B_d,C_x2,C_xinf)
    
    Tilde_Ax = Q*Lambda*Q'/2;
    Tilde_Bu = Q*L_j'*B_u/2;
    Tilde_Bd = Q*L_j'*B_d/2;
    
    if ~isempty(C_x2)
        Tilde_Cy = C_x2*L_j*inv(L_j'*L_j)*Q'/2;
    else
        Tilde_Cy=[];
    end
    if ~isempty(C_xinf)
        Tilde_Cz = C_xinf*L_j*inv(L_j'*L_j)*Q'/2;
    else
        Tilde_Cz=[];
    end
end

