function [optimal_Kpf,optimal_WD,optimal_XD] =...
    Optimal_Peva_D_stability_H2_Hinf_control(...
        B_u,B_d,C_x2,D_u2,C_xinf,D_uinf,D_dinf,...
        Lambda,L_j,Q,...
        alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P,...
        c_H2,c_Hinf,Yalmip_sdpsettings)
%   A) Objective:
%
%   Function to perform partial allocation of eigenvalues with the goal of 
%   optimal H2 and/or H∞ control, ensuring D-stability delimited by 
%   vertical-strip, sector, and disc regions.
%
%   B) Requirements
%   
%       Must have data packages installed:
%
%       1) Control System Toolbox
%       2) Robust Control Toolbox
%       3) YALMIP - MATLAB toolbox for optimization modeling 
%
%       Matrices "B_u","B_d","C_x2","D_u2","C_xinf","D_uinf","D_dinf" must 
%           be real
%
%       Scalar numbers "alpha_v","beta_v","tetha_s","r_d","q_d","c_H2",
%           "c_Hinf" must be positive real
%
%   C) Inputs:
%
%       1) "B_u","B_d","C_x2","D_u2","C_xinf","D_uinf","D_dinf" -> Matrices  
%           that constitute the original state space model, simplified as 
%           follows:
%
%             d(x(t)) = A_x*x(t)+B_u*u(t)+B_d*d(t)
%                y(t) = C_x2*x(t)+D_u2*u(t)
%                z(t) = C_xinf*x(t)+D_uinf*u(t)+D_dinf*d(t)
%
%               x(t) -> state-space vector.
%            d(x(t)) -> derivate of state-space vector.
%               u(t) -> control-input vector.
%               d(t) -> disturb-input vector.
%               y(t) -> vector of undisturbed outputs.
%               z(t) -> vector of disturbed outputs.
%
%       2) "Lambda"-> the diagonal matrix of eigenvalues
%   
%       3) "L_j" -> Matrix  whose columns are left eigenvectors linearly 
%                   independent. 
%       
%       4) "Q" transformation matrix.
%
%       5) "alpha_v" and "beta_v" -> Maximum and minimum real real values
%               that delimite D-stable Vertical-Strip region. if the data 
%               is an empty array, then there will be no allocation to 
%               D-stable region of type X.
%
%       6) "tetha_s" -> Angle of D-stable Sector Region. if the data 
%               is an empty array, then there will be no allocation to 
%               D-stable region of type X.
%
%       7) "r_d","q_d" -> Radius and Center of D-stable Disk Region. if the 
%               data is an empty array, then there will be no allocation to 
%               D-stable region of type X.
%
%       8) "w_H" -> Maximum absolute value for the complex part of a 
%                   complex pole pair:
%                   abs{imag{s}}<=w_H
%
%       9) "e_P" -> Damping parameter of the stability parabola:
%                   imag{s}^2<= -e_P*real{s}
%
%       8) "cH2", "cHinf" -> weights of the H2 norms of the cost function 
%               used in the optimization.if the data is an empty array, 
%               then there will be no optimization of the associated norm.
%
%       9) "Yalmip_sdpsettings" -> data to configure the options related 
%               to the "sdpsettings" command present in "YALMIP". Please 
%               visit the "YALMIP" website for more details on how to 
%               configure "sdpsettings".
%
%   D) Outputs:
%
%       Reduced system described according to equations:
%
%        L_j'*x(t) = h{x}(t) 
%       d(h{x}(t)) = Q*Lambda*Q'*h{x}(t)+...
%                    Q*L_j'*B_u*u(t)+...       
%                    Q*L_j'*B_d*d(t)
%             y(t) = C_x2*L_j*inv(L_j'*L_j)*Q'*h{x}(t)+...
%                    D_u2*u(t)
%             z(t) = C_xinf*L_j*inv(L_j'*L_j)*Q'*h{x}(t)+...
%                    D_uinf*u(t)+D_dinf*d(t)
%             u(t) = Kpf*x(t)
%              Kpf = KD*L_j'*Q 
%
%       "optimal_WD" -> Value of Real Matrix WD is solution of SDP 
%                       optimization
%       "optimal_XD" -> Value of Symmetric Matrix XD is solution 
%                       of SDP optimization. 
%       "optimal_YH2" -> Value of Symmetric Matrix YH2 is solution 
%                        of SDP optimization.
%       "optimal_pH2" -> Value of Symmetric Matrix pH2 is solution 
%                        of SDP optimization.
%       "optimal_pHinf" -> Value of Symmetric Matrix pHinf is solution 
%                          of SDP optimization.
%--------------------------------------------------------------------------
%%  Compute Constants 
%
    n_u = size(B_u,2);
    n_d = size(B_d,2);
    n_y = size(C_x2,1);
    n_z = size(C_xinf);
    n_ola = size(Lambda,1);
    
    Tilde_Lambda = Q*Lambda*Q';
    
    Tilde_Bu = Q*L_j'*B_u;
    
    Tilde_Bd = Q*L_j'*B_d;
    
%%  Make LMIs

    XD  = sdpvar(n_ola,n_ola,'symmetric');
    WD  = spdvar(n_u,n_ola,'full');
    LMIs = XD>=eps*eye(nola);
    
    if ~isempty(alpha_v)
        Tilde_Lambda_Valpha = Lambda+alpha_v*eye(n_ola);
        LMIs = [ LMIs, Tilde_Lambda_Valpha*XD+XD*Tilde_Lambda_Valpha'+...
                 Tilde_Bu*WD+WD'*Tilde_Bu<=-eps*eye(n_ola)];
    end
    
    if ~isempty(beta_v)
        Tilde_Lambda_Dqr = Lambda+q_d*eye(n_ola); 
        LMIs = [ LMIs, Tilde_Lambda_Vbeta*XD+XD*Tilde_Lambda_Vbeta'+...
                 Tilde_Bu*WD+WD'*Tilde_Bu<=-eps*eye(n_ola) ];
    end
    
    if ~isempty(r_d)
        LMIs = [ LMIs,...
                [ -r_d*XD, Tilde_Lambda_Dqr*XD+Tilde_Bu*WD;...
                  XD*Tilde_Lambda_Dqr'+WD'*Tilde_Bu', -r_d*XD ]...
                  <=eps*eye(2*n_ola)...
               ];
    end
    
    if ~isempty(theta_s)
        Tilde_Lambda_sin = Lambda*sin(theta_s);
        Tilde_Lambda_cos = Lambda*cos(theta_s);
        Tilde_Bu_sin = Tilde_Bu*sin(theta_s);
        Tilde_Bu_cos = Tilde_Bu*cos(theta_s);
        LMIs =[ LMIs,...
                [ Tilde_Lambda_sin*XD+XD*Tilde_Lambda_sin'+...
                  Tilde_Bu_sin*WD+WD'*Tilde_Bu_sin',...
                  Tilde_Lambda_cos*XD-XD*Tilde_Lambda_cos'+...
                  Tilde_Bu_cos*WD-WD'*Tilde_Bu_cos';...
                  XD*Tilde_Lambda_cos'-Tilde_Lambda_cos*XD+...
                  WD'*Tilde_Bu_cos'-Tilde_Bu_cos*WD,...
                  Tilde_Lambda_sin*XD+XD*Tilde_Lambda_sin'+...
                  Tilde_Bu_sin*WD+WD'*Tilde_Bu_sin'...
                ]<=eps*eye(2*n_ola)...
           ];
    end
    
    if ~isempty(w_H)
        LMIs = [ LMIs,...
                   [   -w_H*XD,...
                        XD*Tilde_Lambda'-Tilde_Lambda*XD+...
                        WD'*Tilde_Bu'-Tilde_Bu*WD;...
                        Tilde_Lambda*XD-XD*Tilde_Lambda'+...
                        Tilde_Bu*WD-WD'*Tilde_Bu',...
                        -w_H*XD...
                    ]<=eps*eye(2*n_ola)...
                ]; 
    end
    
    if ~isempty(e_P)
        Tilde_Lambda_Pe = (Lambda-e_P*eye(n_ola))/2;                
        LMIs = [ LMIs,...
                    [   Tilde_Lambda_Pe*XD+XD*Tilde_Lambda_Pe'+...
                        Tilde_Bu/2*WD+WD'*Tilde_Bu'/2,...
                        -Tilde_Lambda*XD-Tilde_Bu*WD;...
                        -XD*Tilde_Lambda'-WD'*Tilde_Bu',...
                        Tilde_Lambda*XD/2+XD*Tilde_Lambda'/2+...
                        Tilde_Bu/2*WD+WD'*Tilde_Bu'/2 ...
                    ]<=eps*eye(2*n_ola)...
                ];    
    end
    
    if ~isempty(c_H2)
        
        Tilde_Cx2 = C_x2*L_j/(L_j'*L_j);
        Tilde_Cx2 = Tilde_Cx2*Q;
        YH2 = sdpvar(n_y,n_ola,'full');
        pH2 = sdpvar(1,1,'symmetric');   
        LMIs = [ LMIs, Tilde_Lambda*XD+XD*Tilde_Lambda'+...
                 Tilde_Bu*WD+WD'*Tilde_Bu+...
                 Tilde_Bd*Tilde_Bu'<=-eps*eye(n_ola) ];
        LMIs = [ LMIs,...
                    YH2>=eps*eye(n_y),...
                    trace(YH2)<=pH2,...
                    pH2>=eps,...
                    [   -YH2, Tilde_Cx2*XD+D_u2*WD;...
                         XD*Tilde_Cx2'+WD'*D_u2, -XD...
                    ]<=eps*eye(n_y+n_x)...
                ];
    end
    
    if ~isempty(c_Hinf)
        Tilde_Cxinf = C_xinf*L_j/(L_j'*L_j);
        Tilde_Cxinf = Tilde_Cxinf*Q;
        pHinf = sdpvar(1,1,'symmetric');   
        LMIs = [ LMIs,...
                    pHinf>=eps,...
                    [   Tilde_Lambda*XD+XD*Tilde_Lambda'+...
                        Tilde_Bu*WD+WD'*Tilde_Bu,...
                        Tilde_Bd,...
                        XD*Tilde_Cxinf'+WD*D_uinf';...
                        Tilde_Bd', -pHinf*eye(n_d), D_dinf';...
                        Tilde_Cxinf*XD+D_uinf*WD,...
                        D_dinf, -pHinf*eye(n_z)...
                    ]<=eps*eye(n_x+n_d+n_z)...
                ];
    end
%% Do optimization

    if ~isempty(c_H2)
         if ~isempty(c_Hinf)
            optimize(LMI_first_step,...
                     c_H2*pH2+c_Hinf*pHinf,...
                     Yalmip_sdpsettings); 
         else
             optimize(LMI_first_step,...
                      c_H2*pH2,...
                      Yalmip_sdpsettings);
         end
    else
        if ~isempty(c_Hinf)
            optimize(LMI_first_step,...
                     c_Hinf*pHinf,...
                     Yalmip_sdpsettings); 
         else
             optimize(LMI_first_step,...
                      [],...
                      Yalmip_sdpsettings);
         end
    end  

%% Compute Feedback Matrix

    optimal_WD = value(WD);
    optimal_XD = value(XD);
    optimal_Kpf=optimal_WD/optimal_XD;
    optimal_Kpf=optimal_Kpf*L_j'*Q;
end

