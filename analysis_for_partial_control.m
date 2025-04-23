function [rank_Obsv,rank_Ctrb,YTX,Lambda_1] =...
    analysis_for_partial_control(A_x,B_u,B_d,...
        C_x2,C_xinf,D_u2,D_uinf,D_dinf,...
        Left_Eigenvectors,Eigenvalues,Right_Eigenvectors)
    Lambda_1=Left_Eigenvectors'*A_x*Right_Eigenvectors;
    YTX=Left_Eigenvectors'*Right_Eigenvectors;
    rank_Obsv = rank(Lambda_1);
    rank_Ctrb = inputArg2;
end

