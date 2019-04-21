function [Net_a,Net_f] = Anova(genes, regulators, expressiondata)

    warning ('off','all')

    nr = size(regulators,1);
    tfs = zeros(1,nr);

    parfor i = 1:nr
        try
            tfs(1,i) = find(strcmp(genes,regulators{i}));
        catch
            continue
        end
    end

    tfs = tfs(tfs ~= 0);

    ngenes = size(expressiondata,1);
    ntf = size(tfs,2);

    n2_a=zeros(ngenes,ntf);
    n2_f=zeros(ngenes,ntf); 

    tfexpression = expressiondata(tfs,:);

    for i = 1:ngenes

        expredatap_2 = expressiondata(i,:);

        parfor j = 1:ntf   

            expredatap_1 = tfexpression(j,:); 

            [p_a,SS_a] = anova2([expredatap_1;expredatap_2],1,'off');
            [p_f,SS_f] = friedman([expredatap_1;expredatap_2],1,'off');

            if p_a(1,2) < 0.05 
                np_a=SS_a{3,2}/SS_a{5,2};   %SSA/SStotal
            else
                np_a=0;
            end
            if p_f < 0.05 
                np_f=SS_f{2,2}/SS_f{4,2};   %SSA/SStotal
            else 
                np_f=0;
            end

            expredatan_1 = -expredatap_1;
            expredatan_2 = expredatap_2;

            [p_a,SS_a] = anova2([expredatan_1;expredatan_2],1,'off');
            [p_f,SS_f] = friedman([expredatan_1;expredatan_2],1,'off')

            if p_a(1,2) < 0.05 
                nn_a=SS_a{3,2}/SS_a{5,2};  %SSA/SStotal
            else
                nn_a=0;
            end
            if p_f < 0.05 
                nn_f=SS_f{2,2}/SS_f{4,2};  %SSA/SStotal
            else 
                nn_f=0;
            end

            if np_a>nn_a
                n2_a(i,j)=np_a;
            else
                n2_a(i,j)=-nn_a;
            end

            if np_f>nn_f
                n2_f(i,j)=np_f;
            else
                n2_f(i,j)=-nn_f;
            end
        end
    end

    Net_a = cell(nnz(n2_a),3);
    
    r=1;

    for i = 1:ntf
        for j = 1:ngenes 
            if n2_a(j,i) == 0 
                continue
            else
                Net_a{r,1} = genes{tfs(i)};
                Net_a{r,2} = genes{j};
                Net_a{r,3} = n2_a(j,i);
                r = r+1;
            end
        end
    end
    
    Net_f = cell(nnz(n2_f),3);
    
    r=1;

    for i = 1:ntf
        for j = 1:ngenes 
            if n2_f(j,i) == 0 
                continue
            else
                Net_f{r,1} = genes{tfs(i)};
                Net_f{r,2} = genes{j};
                Net_f{r,3} = n2_f(j,i);
                r = r+1;
            end
        end
    end
    
end


