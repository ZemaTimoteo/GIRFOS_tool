function bound = sharp_imProfile(auxTest,fig)
    

    % division in two parts
    part1_vect = auxTest(1:round(size(auxTest,1)/2),:);
    part2_vect = auxTest(round(size(auxTest,1)/2)+1:end,:);
    
    % figures
    if fig == 'True'
        figure()
        subplot(121)
        plot(auxTest)
        subplot(122)
        plot(part1_vect)
        hold on
        plot(part2_vect)
    end
    
    % --- part1 analysis ---
    upLIM_p1   = 0.8*max(part1_vect);
    aux_UP_p1  = find(part1_vect>=upLIM_p1);
    difUP_p1    = diff(aux_UP_p1);  % account for min relative
    auxDifUP_p1 = find(difUP_p1>1); % account for min relative
    if isempty(auxDifUP_p1) 
        auxDifUP_p1 = 0;
    end
        
    lowLIM_p1    = 0.2*max(part1_vect);
    aux_LOW_p1   = find(part1_vect<=lowLIM_p1);
    difLOW_p1    = diff(aux_LOW_p1);  % account for min relative
    auxDifLOW_p1 = find(difLOW_p1>1); % account for min relative
    if isempty(auxDifLOW_p1)
        auxDifLOW_p1 = size(aux_LOW_p1,1);
    end

    firstUP    = aux_UP_p1(auxDifUP_p1(end)+1);
    firstLOW   = aux_LOW_p1(auxDifLOW_p1(end));
    
    % --- part2 analysis ---
    upLIM_p2    = 0.8*max(part2_vect);
    aux_UP_p2   = find(part2_vect>=upLIM_p2);
    difUP_p2    = diff(aux_UP_p2);  % account for min relative
    auxDifUP_p2 = find(difUP_p2>1); % account for min relative
    if isempty(auxDifUP_p2) 
        auxDifUP_p2 = 0;
    end
    
    lowLIM_p2    = 0.2*max(part2_vect);    
    aux_LOW_p2   = find(part2_vect<=lowLIM_p2);
    difLOW_p2    = diff(aux_LOW_p2);  % account for min relative
    auxDifLOW_p2 = find(difLOW_p2>1); % account for min relative
    if isempty(auxDifLOW_p2)
        auxDifLOW_p2 = 0;
    end   
        
    secondUP    = aux_UP_p2(end-auxDifUP_p2(1));
    secondLOW   = aux_LOW_p2(auxDifLOW_p2(end)+1);
    
    

    % return output
    bound.firstUP = firstUP;
    bound.secondUP = secondUP;
    bound.firstLOW = firstLOW;
    bound.secondLOW = secondLOW;

end
