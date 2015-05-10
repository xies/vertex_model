function tisArr = assemble_model(DIR,n)

    tisArr(1:n) = Tissue;
    
    for i = 1:n
        load([DIR '/model_step_' num2str(i) '.mat'])
        tisArr(i) = tis;
    end

end