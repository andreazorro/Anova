function [sNet_f] = Friedman(genes, regulators, expressiondata)

    warning ('off','all')

    expressiondata = expressiondata';
    
    tfs = cellfun(@(x)find(strcmp(x,gene_names)),regulators,'UniformOutput',false);
    tfs(logical(cellfun('isempty',tfs))) = {0};
    tfs = cell2mat(tfs);
    tfs = tfs(tfs ~= 0);

    ngenes = size(expressiondata,1);
    ntf = size(tfs,1);

    n2_f=zeros(ngenes,ntf); 

    tfexpression = expressiondata(tfs,:);

    for i = 1:ngenes

        expredatap_2 = expressiondata(i,:);

        parfor j = 1:ntf   
            
            if i==tfs(1,j)
                continue
            else
                expredatap_1 = tfexpression(j,:); 

                [p_f,SS_f] = friedman([expredatap_1;expredatap_2],1,'off');

                if p_f < 0.05 
                    np_f=SS_f{2,2}/SS_f{4,2};   %SSA/SStotal
                else 
                    np_f=0;
                end

                expredatan_1 = -expredatap_1;
                expredatan_2 = expredatap_2;

                [p_f,SS_f] = friedman([expredatan_1;expredatan_2],1,'off')

                if p_f < 0.05 
                    nn_f=SS_f{2,2}/SS_f{4,2};  %SSA/SStotal
                else 
                    nn_f=0;
                end

                if np_f>nn_f
                    n2_f(i,j)=np_f;
                else
                    n2_f(i,j)=-nn_f;
                end
            end
        end
    end

    n2_f=abs(n2_f);   
   
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

    sNet_f = sortrows(Net_f,3,'descend');
    
end


