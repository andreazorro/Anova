function [sNet_a] = Anova(genes, regulators, expressiondata)

    warning ('off','all')

    expressiondata = expressiondata';
    
    tfs = cellfun(@(x)find(strcmp(x,gene_names)),regulators,'UniformOutput',false);
    tfs(logical(cellfun('isempty',tfs))) = {0};
    tfs = cell2mat(tfs);
    tfs = tfs(tfs ~= 0);

    ngenes = size(expressiondata,1);
    ntf = size(tfs,1);

    n2_a=zeros(ngenes,ntf);
  
    tfexpression = expressiondata(tfs,:);

    for i = 1:ngenes

        expredatap_2 = expressiondata(i,:);

        parfor j = 1:ntf   
            
            if i==tfs(1,j)
                continue
            else
                expredatap_1 = tfexpression(j,:); 

                [p_a,SS_a] = anova2([expredatap_1;expredatap_2],1,'off');
                
                if p_a(1,2) < 0.05 
                    np_a=SS_a{3,2}/SS_a{5,2};   %SSA/SStotal
                else
                    np_a=0;
                end
               
                expredatan_1 = -expredatap_1;
                expredatan_2 = expredatap_2;

                [p_a,SS_a] = anova2([expredatan_1;expredatan_2],1,'off');
         
                if p_a(1,2) < 0.05 
                    nn_a=SS_a{3,2}/SS_a{5,2};  %SSA/SStotal
                else
                    nn_a=0;
                end
               
                if np_a>nn_a
                    n2_a(i,j)=np_a;
                else
                    n2_a(i,j)=-nn_a;
                end
            end
        end
    end

    n2_a=abs(n2_a);
   
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

    sNet_a = sortrows(Net_a,3,'descend');
      
end


