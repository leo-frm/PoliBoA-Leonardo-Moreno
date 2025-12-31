function doGet() {
  try {
    // Conecta à planilha "Dados Site"
    const spreadsheet = SpreadsheetApp.openById('1e7fIThPQbYtwPf8kN81l0faz3obTY8Ouk0ULbcKpwcw');
    const sheet = spreadsheet.getSheetByName("dadosSite");
    
    if (!sheet) {
      throw new Error("A aba 'dadosSite' não foi encontrada na planilha.");
    }
    
    const allData = sheet.getDataRange().getValues();
    const headers = allData.shift(); // Remove cabeçalhos
    
    // Obtém hora atual no timezone de São Paulo
    const now = new Date();
    const currentDate = Utilities.formatDate(now, "America/Sao_Paulo", "dd/MM/yyyy");
    const currentTime = Utilities.formatDate(now, "America/Sao_Paulo", "HH:mm:ss");
    const currentHour = parseInt(Utilities.formatDate(now, "America/Sao_Paulo", "HH"), 10);
    
    console.log("🔍 Data atual (SP):", currentDate);
    console.log("🕐 Hora atual (SP):", currentTime, "(" + currentHour + "h)");
    console.log("📊 Total de linhas:", allData.length);

    // Índices das colunas (conforme estrutura fornecida)
    const dataIdx = headers.indexOf('data');
    const horaIdx = headers.indexOf('hora');
    
    let latestDataRow = null;
    let latestHour = -1;
    let latestDate = null;
    let foundTodayData = false;

    // Procura dados de HOJE que já passaram
    for (let i = allData.length - 1; i >= 0; i--) {
      const row = allData[i];
      const rowDate = row[dataIdx];
      const rowHourStr = row[horaIdx];
      
      if (!rowDate || !rowHourStr) continue;
      
      // Normaliza a data para dd/MM/yyyy
      let rowDateStr = '';
      if (rowDate instanceof Date) {
        rowDateStr = Utilities.formatDate(rowDate, "America/Sao_Paulo", "dd/MM/yyyy");
      } else {
        rowDateStr = String(rowDate).trim();
      }
      
      // Extrai a hora
      let rowHour = -1;
      if (rowHourStr instanceof Date) {
        rowHour = rowHourStr.getHours();
      } else {
        const hourStr = String(rowHourStr).trim();
        const hourMatch = hourStr.match(/^(\d{1,2})/);
        if (hourMatch) {
          rowHour = parseInt(hourMatch[1], 10);
        }
      }
      
      if (isNaN(rowHour)) continue;
      
      // Debug das últimas 10 linhas
      if (i >= allData.length - 10) {
        console.log(`Linha ${i+2}: Data="${rowDateStr}" Hora=${rowHour}h`);
      }
      
      // Prioriza dados de HOJE que já passaram
      if (rowDateStr === currentDate && rowHour <= currentHour) {
        if (rowHour > latestHour) {
          latestDataRow = row;
          latestHour = rowHour;
          latestDate = rowDateStr;
          foundTodayData = true;
          console.log(`✅ Candidato atualizado: ${latestDate} ${String(latestHour).padStart(2, '0')}:00`);
        }
      }
    }
    
    if (foundTodayData) {
      console.log(`✅ DADOS FINAIS DE HOJE: ${latestDate} ${String(latestHour).padStart(2, '0')}:00`);
    }

    // Se não achou dados de hoje, pega o registro mais recente disponível
    if (!foundTodayData) {
      console.log(`⚠️ Sem dados para hoje até as ${currentHour}h. Buscando registro mais recente...`);
      
      // Pega simplesmente a última linha válida
      for (let i = allData.length - 1; i >= 0; i--) {
        const row = allData[i];
        const rowDate = row[dataIdx];
        const rowHourStr = row[horaIdx];
        
        if (!rowDate || !rowHourStr) continue;
        
        let rowDateStr = '';
        if (rowDate instanceof Date) {
          rowDateStr = Utilities.formatDate(rowDate, "America/Sao_Paulo", "dd/MM/yyyy");
        } else {
          rowDateStr = String(rowDate).trim();
        }
        
        let rowHour = -1;
        if (rowHourStr instanceof Date) {
          rowHour = rowHourStr.getHours();
        } else {
          const hourMatch = String(rowHourStr).trim().match(/^(\d{1,2})/);
          if (hourMatch) {
            rowHour = parseInt(hourMatch[1], 10);
          }
        }
        
        if (isNaN(rowHour)) continue;
        
        latestDataRow = row;
        latestHour = rowHour;
        latestDate = rowDateStr;
        break; // Pega a primeira linha válida de trás para frente
      }
      
      if (latestDate) {
        console.log(`✅ Usando dados de: ${latestDate} ${String(latestHour).padStart(2, '0')}:00`);
      }
    }

    if (!latestDataRow) {
      throw new Error(`Nenhum dado disponível na planilha`);
    }

    // Função auxiliar para converter valores com vírgula para número
    function parseValue(val) {
      if (val === null || val === undefined || val === '' || val === 'X') return 0;
      if (typeof val === 'number') return val;
      const strVal = String(val).replace(',', '.');
      const parsed = parseFloat(strVal);
      return isNaN(parsed) ? 0 : parsed;
    }

    // Monta resposta com os dados da linha correta
    const responseData = {
      nivelRio: parseValue(latestDataRow[headers.indexOf('nivelRio')]),
      precipitacaoAcumulada: parseValue(latestDataRow[headers.indexOf('precipitacaoAcumulada')]),
      umidadeAr: parseValue(latestDataRow[headers.indexOf('umidadeAr')]),
      temperatura: parseValue(latestDataRow[headers.indexOf('Temperatura')]),
      precipitacaoMM: parseValue(latestDataRow[headers.indexOf('PrecipitacaoMM')]),
      condicaoTempo: latestDataRow[headers.indexOf('condicaoTempo')] || '',
      chanceChuva: parseValue(latestDataRow[headers.indexOf('chanceChuva')]),
      hasPrecipitation: Boolean(latestDataRow[headers.indexOf('hasPrecipitation')]),
      riscoCalculado: latestDataRow[headers.indexOf('riscoCalculado')] || 'Indeterminado',
      dataColeta: `${latestDate} ${String(latestHour).padStart(2, '0')}:00`,
      horaAtualServidor: `${currentDate} ${currentTime}`
    };
    
    console.log("✅ DADOS ENVIADOS:", JSON.stringify(responseData));
    
    return ContentService.createTextOutput(JSON.stringify(responseData))
      .setMimeType(ContentService.MimeType.JSON);

  } catch (err) {
    console.error("❌ Erro na função doGet:", err.message);
    return ContentService.createTextOutput(JSON.stringify({ 
      error: err.message,
      timestamp: new Date().toISOString() 
    }))
    .setMimeType(ContentService.MimeType.JSON);
  }
}

/**
 * Calcula o risco baseado em dados meteorológicos
 * @param {number} precipitacao Precipitação acumulada
 * @param {number} umidade Umidade do ar (%)
 * @param {string} condicao Condição do tempo
 * @param {number} chanceChuva Chance de chuva (%)
 * @param {number} temperatura Temperatura (°C)
 * @param {number} nivelRio Nível do rio
 * @return {string} Nível de risco (Baixo, Médio, Alto ou Crítico)
 * @customfunction
 */
function CALCULAR_RISCO(precipitacao, umidade, condicao, chanceChuva, temperatura, nivelRio) {
  // Validação de entrada - aceita zero como valor válido
  if (precipitacao == null || umidade == null || !condicao || 
      chanceChuva == null || temperatura == null || nivelRio == null) {
    return "Dados incompletos";
  }
  
  // Mapeamento de condições
  const condicoes = {
    "claro": 0.0,
    "ensolarado": 0.0,
    "parcialmente nublado": 0.0,
    "nublado": 0.05,
    "chuvisco leve": 0.1,
    "chuva irregular": 0.2,
    "chuva leve": 0.3,
    "chuva": 0.6,
    "tempestade": 1.0
  };
  
  // Normalizar condição (remover espaços extras e converter para minúscula)
  const condicaoNormalizada = condicao.toString().trim().toLowerCase();
  const condicaoVal = condicoes[condicaoNormalizada] !== undefined ? condicoes[condicaoNormalizada] : 0.4;
  
  // Normalizações
  const precipitacaoNorm = precipitacao / 77.8;
  const umidadeNorm = umidade / 100;
  const chanceChuvaVal = typeof chanceChuva === 'string' && chanceChuva.includes('Sim') ? 100 : chanceChuva;
  const chanceChuvaaNorm = chanceChuvaVal / 100;
  const temperaturaNorm = (temperatura - 5) / (40 - 5);
  const nivelRioNorm = (nivelRio - 100) / (500 - 100);
  
  // Cálculo do risco
  const x = (
    precipitacaoNorm * 0.30 +
    umidadeNorm * 0.025 +
    condicaoVal * 0.075 +
    chanceChuvaaNorm * 0.15 +
    temperaturaNorm * 0.05 +
    (1 - nivelRioNorm) * 0.4
  ) * 100;
  
  // Classificação do risco
  if (x < 25) {
    return "Baixo";
  } else if (x < 50) {
    return "Médio";
  } else if (x < 70) {
    return "Alto";
  } else {
    return "Crítico";
  }
}