fileName = 'error000_modified.dat';

plotFlag = 0;
for errorType = 1:2
    rawData = importdata(fileName);
    rawNumbers = csvread(fileName,1,4);
    
    iterNum = str2double(rawData.textdata(:,1));
    
    order = str2double(rawData.textdata(:,2));
    
    mesh = rawData.textdata(:,4);
    
    nRuns = length(order);
    
    tetNum = zeros(nRuns,1);
    triNum = zeros(nRuns,1);
    
    errorU = zeros(nRuns,1);
    errorDiffU = zeros(nRuns,1);
    
    for i = 2:nRuns
        getTetNum = textscan(mesh{i},'tet_%f_%f_%f');
        if ~isempty(getTetNum{1})
            tetNum(i) = getTetNum{1};
            errorU(i) = rawNumbers(i-1,11);
            errorDiffU(i) = rawNumbers(i-1,16);
            
        end
        getTriNum = textscan(mesh{i},'tri_%f_%f');
        if ~isempty(getTriNum{1})
            triNum(i) = getTriNum{1};
            errorU(i) = rawNumbers(i-1,10);
            errorDiffU(i) = rawNumbers(i-1,14);
        end
        
    end
    
    is1 = order == 1;
    is2 = order == 2;
    is3 = order == 3;
    is4 = order == 4;
    is5 = order == 5;
    is6 = order == 6;
    
    isTet = tetNum ~= 0;
    isTri = triNum ~= 0;
    
    
    isLastIterTet = iterNum == 1e4;
    isLastIterTri = iterNum == 1e6;
    
    isLastIter = isLastIterTet + isLastIterTri;
    
    elem = [tetNum triNum];
    isElem = [isTet isTri];
    isOrder = [is1 is2 is3 is4 is5 is6];
    error = [errorU errorDiffU];
    
    %% Store results in a cell
    
    % Get array of different meshes used
    meshVals(1) = {unique(tetNum(1:end-1))};
    meshVals(2) = {unique(triNum)};
    
    % Get array of different orders used
    orderVals = unique(order);
    
    % Pre-allocate memory in cell
    errorMat = cell(...
        max(numel(meshVals{1}),numel(meshVals{2}))-2,numel(orderVals)-1,2);
    size(errorMat)
    ordMat = cell(max(numel(meshVals{1}),numel(meshVals{2}))-1,...
        numel(orderVals)-1,2);
    
    
    %%
    figureNum = 1;
    
    for eleType = 1:2
        for ordNum = 1:5
            
            %         eleType = 2; % tet = 1; tri = 2;
            %         ordNum = 2;
            %errorType = 2; % in U = 1; in diffU = 2;
            
            elemToPlot = (isOrder(:,ordNum).*isElem(:,eleType)...
                .*isLastIter) ~= 0;
            errorToPlot = error(elemToPlot,errorType);
            runToPlot = elem(elemToPlot,eleType);
            
            if eleType == 2 && ordNum ==5
                errorToPlot = errorToPlot(1:3);
                runToPlot = runToPlot(1:3);
            end
            
            for k = 1:length(runToPlot)
                meshIn = find(runToPlot(k) == meshVals{eleType});
                if ~isempty(meshIn)
                    errorMat{k,ordNum,eleType} = errorToPlot(k);
                end
            end
            
            numMeshes = numel(unique(elem(elemToPlot,eleType)));
            
            % Fit the data
            p = polyfit(log(runToPlot),log(errorToPlot),1);
            
            ordMat{end,ordNum,eleType} = sprintf('%.2f',-p(1))
            
            fitValX = linspace(min(runToPlot),max(runToPlot),10);
            fitValY = exp(p(2))*fitValX.^(p(1));
            
            
            
            if plotFlag
                titleName = ...
                    sprintf(['eleType: %i\nNum runs: %i \nNum meshes: %i\n'...
                    'Order: %i\nReal order: %.3f'], ...
                    eleType,sum(elemToPlot),numMeshes,ordNum,p(1))
                
                figure(figureNum)
                loglog(runToPlot,errorToPlot,'.',fitValX,fitValY,'-')
                title(titleName)
            end
            figureNum = figureNum + 1;
        end
    end
    
    
    %% Get order of accuracy
    
    for elType = 1:2
        for ord = 1:size(errorMat,2)
            for msh = 1:size(errorMat,1)-1
                
                y2 = errorMat{msh+1,ord,elType};
                y1 = errorMat{msh,ord,elType};
                
                x2 = meshVals{elType}(msh+2);
                x1 = meshVals{elType}(msh+1);
                
                if ~isempty(y2) && ~isempty(y1)
                    ordMat{msh+1,ord,elType} = ...
                        sprintf('%.2f',-log(y2/y1)/log(x2/x1));
                end
                
            end
        end
    end
    
    %% Combine error and order of accuracy tables
    for eT = 1:2 % error type
        resMatRows = max(size(errorMat(:,:,eT),1),size(ordMat(:,:,eT),1));
        resMatCols = size(errorMat(:,:,eT),2) + size(ordMat(:,:,eT),2);
        resMat = cell(resMatRows,resMatCols)
        
        for i =1:resMatCols/2
            for j = 1:resMatRows
                try  resMat(j,2*i-1) = errorMat(j,i,eT); end
                try resMat(j,2*i) = ordMat(j,i,eT); end
            end
        end
        
        
        
        
        resMat = resMat';
        
        rowFormat = repmat([ {'\\hline'} {''}],1,resMatRows)
        
        
        
        
        if eT == 1
            nRepeats = 3;
            nMeshes = length(meshVals{eT}(3:end))
        elseif eT == 2
            nRepeats = 2;
            nMeshes = length(meshVals{eT}(2:end))
        end
        
        firstRow = cell(1,nMeshes+1);
        
        for i = 1:nMeshes
            val = repmat(meshVals{eT}(i+1),1,nRepeats);
            firstRow{i} = sprintf([repmat('%ix',1,nRepeats-1) '%i'],val);
        end
        
        firstRow{end} = 'Overall Order of Accuracy';
        firstRow
        
        resMat = [firstRow; resMat]
        
        %%
        if eT == 1
            resMat = [resMat(:,1:4) resMat(:,end)];
        end
        
        % Second column
        secondCol = [{'Mesh:'};repmat([{'$L_2$ error'};{'$\\mathcal{O}(L_2)$'}],resMatCols/2,1)]
        
        
        % First column
        firstCol = [{'Polynomial Order'}];
        
        for i = 1:resMatCols/2
            firstCol = [firstCol; {['\\multirow{2}{*}{$p = '...
                num2str(i) '$}']};{''}];
        end
        
        %;repmat([{'\\multirow{2}{*}{Multirow}'};{''}],resMatCols/2,1)]
        
        %%
        resMat = [firstCol secondCol resMat]
        %%
        if eT == 1
            if errorType == 1
                captionString = ['\\caption{Accuracy of HiFiLES for NS'...
                    ' equations with source term in tetrahedral meshes at $t = 10$.'...
                    ' $L_2$ error is the $L_2$-norm of the error in the energy field: $\\rho e$}\n'...
                    '\\label{table:tetsError' num2str(errorType) '}'];
            elseif errorType == 2
                captionString = ['\\caption{Accuracy of HiFiLES for NS'...
                    ' equations with source term in tetrahedral meshes at $t = 10$.'...
                    ' $L_2$ error is the $L_2$-norm of the error in the gradient of the energy field:'...
                    '$\\frac{\\partial}{\\partial x_i} (\\rho e)$}\n'...
                    '\\label{table:tetsError' num2str(errorType) '}'];
            end
        elseif eT == 2
            if errorType == 1
                captionString = ['\\caption{Accuracy of HiFiLES for NS'...
                    ' equations with source term in triangular meshes at $t = 1$.'...
                    ' $L_2$ error is the $L_2$-norm of the error in the energy field: $\\rho e$}\n'...
                    '\\label{table:trisError' num2str(errorType) '}'];
            elseif errorType == 2
                captionString = ['\\caption{Accuracy of HiFiLES for NS'...
                    ' equations with source term in triangular meshes at $t = 1$.'...
                    ' $L_2$ error is the $L_2$-norm of the error in the gradient of the energy field:'...
                    '$\\frac{\\partial}{\\partial x_i} (\\rho e)$}\n'...
                    '\\label{table:trisError' num2str(errorType) '}'];
            end
        end
        
        tableString=LatexTable(resMat,'%.2e',captionString,...
            repmat({' c'},1,size(resMat,2)),...
            rowFormat);
        
        writeFlag = 1;
        
        tableName = sprintf('summaryTable_ele%i_err%i.tex',eT,errorType);
        
        if writeFlag
            % Write in file
            fileID = fopen(['/home/mlopez14/Desktop' ...
                '/paper/HiFiLES/paper/aiaa_atlanta_2014/' ...
                tableName ],'w');
            fprintf(fileID,tableString);
            fclose(fileID);
            
        end
        
    end
end