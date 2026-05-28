function Off = Insert(Population) 

    len = length(Population);

    Off = [];
    p1 = Population(1);
    for i = 2:len 
        Off = [Off; (p1.decs + Population(i).decs) ./ 2;];
    end 
end