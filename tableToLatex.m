% Função para gerar código LaTeX a partir de uma tabela MATLAB
% Função para gerar código LaTeX a partir de uma tabela MATLAB
function latexTable = tableToLatex(T)
    % Obter os nomes das colunas
    columnNames = T.Properties.VariableNames;
    
    % Criar o formato de alinhamento das colunas, assumindo que todas sejam centradas
    colAlignment = repmat('c', 1, numel(columnNames));
    
    % Inicializar o código LaTeX
    latexTable = ['\\begin{tabular}{|' colAlignment '|} \\hline\n'];
    
    % Adicionar cabeçalho
    latexTable = [latexTable, strjoin(columnNames, ' & '), ' \\\\ \\hline\n'];
    
    % Adicionar dados das linhas
    for i = 1:height(T)
        rowData = T{i, :}; % Obter uma linha da tabela
        latexTable = [latexTable, strjoin(string(rowData), ' & '), ' \\\\ \n'];
    end
    
    % Fechar a tabela
    latexTable = [latexTable, '\\hline\n\\end{tabular}'];
end
