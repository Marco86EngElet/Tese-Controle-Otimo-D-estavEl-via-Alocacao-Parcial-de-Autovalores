%% A) Computar as Normas "H_2" e "H_infinito" do sistema em malha aberta     

    %A.1) Norma H_2 do sistema canonico controlavel de Gody(s)
    
    norma_H2_Godys=norm(Canon_Godys,2);
    
    %A.2) Norma H_infinito do sistema canonico controlavel de Godz(s)
    
    norma_Hinf_Godzs=norm(Canon_Godzs,'inf');
    
%% B) Computar Os polos do sistema em malha aberta 
    
    Polos_Malha_Aberta=cplxpair(pole(Canon_Godys));
