function H = Criar_Funcao_Transferencia_SISO_2nd(K, zeta, omega_n)
    % Criar a funcao de transferencia de segunda ordem
    num = K * omega_n^2; % Numerador
    den = [1, 2*zeta*omega_n, omega_n^2]; % Denominador
    H = tf(num, den); % Criar a funcao de transferencia
end
