function [  x_horizontal_strip,y_inferior_horizontal_strip,...
            y_superior_horizontal_strip,...
            x_sector,y_inferior_sector,y_superior_sector,...
            x_parabola,y_inferior_parabola,y_superior_parabola,...
            x_disk,y_disk,...
            x_alphav,y_alphav,...
            x_betav,y_betav ] = linhas_D_regioes(...
            autovalores,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P)

%A) inicializar vetor para fronteira da LMI Faixa Horizontal
x_horizontal_strip=[];
y_inferior_horizontal_strip=[];
y_superior_horizontal_strip=[];

%B) inicializar vetor para fronteira da LMI setor
x_sector=[] ;
y_inferior_sector=[] ;
y_superior_sector=[] ;

%C) inicializar vetor para fronteira da LMI parabola
x_parabola=[] ;
y_inferior_parabola=[] ;
y_superior_parabola=[] ;

%D) inicializar vetor para fronteira da LMI Disco
x_disk=[] ;
y_disk=[] ;

%E) inicializar vetor para fronteira da LMI Faixa vertical
x_alphav=[] ;
y_alphav=[] ;
x_betav=[] ;
y_betav=[] ;            

%F) Computar a maior e menor parte real dos autovalores do sistema            
min_real_aut=min(real(autovalores)); 
max_real_aut=max(real(autovalores));                

%G) Computar a maior e menor parte real da fronteira LMI Disco
min_real_disk=-q_d-r_d;
max_real_disk=-q_d+r_d;


%H) Computar a menor parte real do conjunto de coordenadas
min_real=min([min_real_aut,min_real_disk,-beta_v]);


%I) Computar a maior parte real do conjunto de coordenadas
max_real=max([max_real_aut,max_real_disk,-alpha_v]);
max_real=min([max_real,0]);


%J) Computar a maior parte imaginaria dos autovalores do sistema
max_imag_aut=max(imag(autovalores));

%K) Computar a maior parte imaginaria da LMI setor
max_imag_sector=-min_real*tan(theta_s);

%M) Computar a maior parte imaginaria da LMI parabola
max_imag_parabola=sqrt(-e_P*min_real);

%N) Computar a maior parte imaginaria do conjunto de coordenadas
max_imag=max([max_imag_aut,max_imag_sector,max_imag_parabola,w_H,r_d]);

%O) Computar todos os cem pontos do vetor para 
for i=1:100
    
    %O.1) Valores no eixo das abscissas
    xt=min_real+(i-1)*(max_real-min_real)/99;

    %O.2) Coordenadas fronteira da LMI Faixa Horizontal
    if ~isempty(w_H)
        x_horizontal_strip(i,1)=xt;
        y_inferior_horizontal_strip(i,1)=-w_H;
        y_superior_horizontal_strip(i,1)=w_H;
    end

    %O.3) Coordenadas fronteira da LMI Setor
    if ~isempty(theta_s)
        x_sector(i,1)=xt;
        y_superior_sector(i,1)=xt*tan(theta_s);
        y_inferior_sector(i,1)=xt*tan(-theta_s);
    end

    %O.4) Coordenadas fronteira da LMI Parabola
    if ~isempty(e_P)
        x_parabola(i,1)=xt;
        y_superior_parabola(i,1)=sqrt(-e_P*xt);
        y_inferior_parabola(i,1)=-sqrt(-e_P*xt);
    end

    %O.5) Coordenadas fronteira da LMI Disco
    if ~isempty(r_d)
        if isempty(q_d)
            x_disk(i,1) = r_d*cos(2*pi*i/100);
            y_disk(i,1) = r_d*sin(2*pi*i/100);
        else
            x_disk(i,1) = -q_d+r_d*cos(2*pi*i/100);
            y_disk(i,1) = -q_d+r_d*sin(2*pi*i/100);
        end
    end

    %O.6) Coordenadas fronteira da LMI Faixa Vertical
    if ~isempty(alpha_v)
        x_alphav(i,1) = -alpha_v;
        y_alphav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
    if ~isempty(beta_v)
        x_betav(i,1) = -beta_v;
        y_betav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
end

end

