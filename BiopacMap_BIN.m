function saida =  BiopacMap_BIN()
global Arquivo Path Sinal Configuracao;

Parametros = [];
Modo = 'Dialogo';
Modo = strcmp(upper(Modo),'DIALOGO');

Resposta = 1;
% PARAMETROS assume o formato { 'Path'  'Arquivo' 'Calibrados'_OU_'Volts' } para as funções AbrirDAS_BIN, AbrirDAS_TXT.

   [Arquivo1 Path1] = uigetfile('*.bin','Importar Arquivo Binário do DAS','MultiSelect','on');


% Atualiza variaveis
%LimparTudo;
Path    = Path1;

for i = 1:length(Arquivo1)
    Arquivo = Arquivo1{i};
    
    %BioMecanica('wait');                             % linha 102 comentada!!!
    FileID=fopen([Path Arquivo],'r','b');
    % Carrega cabecalho da primeira linha - Numero de canais, frequencia de amostragem
    frewind(FileID);
    dim=fread(FileID,2,'uint32');
    temp=fread(FileID,[dim(2) dim(1)],'float32')';
    famost=temp(1,1);
    ncanais=temp(1,2);
    graupoli=dim(2);
    temp=temp(2:end,:);
    % Carrega Linha de base
    BaseLine=temp(1:ncanais,1);
    Configuracao.LinhaBaseSinais = temp(1:ncanais,1);
    Configuracao.RuidoLinhaBaseSinais = temp(1:ncanais,2);
    % Carrega Ganhos positivos
    temp=temp(ncanais+1:end,:);
    PositiveGain=[fliplr(temp(1:ncanais,:)) zeros(ncanais,1)];
    Configuracao.GanhoPositivo = PositiveGain;
    % Carrega Ganhos Negativos
    temp=temp(ncanais+1:end,:);
    NegativeGain=[fliplr(temp(1:ncanais,:)) zeros(ncanais,1)];
    Configuracao.GanhoNegativo = NegativeGain;
    temp=temp(2:end,:);
    temp=temp(ncanais+1:end,:);
    % Carrega Nome de Cada Canal
    dim=fread(FileID,1,'uint32');
    Sinal.Nome = { };
    for j = 1:dim
        Sinal.Nome(j) = {char(fread(FileID,fread(FileID,1,'uint32'),'uchar')')};
    end;
    % Carrega a hora de gravacao
    Hora=char(fread(FileID,fread(FileID,1,'uint32'),'uchar')');
    if isempty(findstr(':',Hora))Linha
        Hora = [];
    end;
    % Carrrega arquivo
    Sinal.Dado=fread(FileID,[ncanais inf],'float32');
    fclose(FileID);
    
     
    % Tira linha de base
    for k = 1:ncanais
        Sinal.Dado(:,k) = (Sinal.Dado(:,k) - BaseLine(k));
    end;
    
    Sinal.Mapa{i} = Sinal.Dado;
end
% Salva configuracao
Configuracao.FrequenciaAmostragem = famost;
Configuracao.ArquivoOriginal = Arquivo;
if ~isempty(Hora)
    Configuracao.Hora = datestr(datenum(Hora),13);
end;

saida.Sinal = Sinal;
saida.Configuracao = Configuracao;
