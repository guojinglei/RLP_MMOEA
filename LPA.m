function [ Labelnew ] = LPA( adjacent_matrix,label )

    N = size(adjacent_matrix,1);
    
    Label1 = label;
    Label2 = Label1;
    Labelnew = Label1;
    flag=1;
    max_iter = 100;
    iter = 0; 
    while(1)
        for i=1:N
            nb_lables = Labelnew(adjacent_matrix(i,:));
            if size(nb_lables,2)>0
                x = tabulate(nb_lables);
                max_nb_labels = x(x(:,2)==max(x(:,2)),1);
                Labelnew(i) = min(max_nb_labels);
            end
        end

        if all(Labelnew==Label1)||all(Labelnew==Label2)
            break;
        else
            if flag==1
                Label1 = Labelnew;
                flag=0;
            else
                Label2 = Labelnew;
                flag=1;
            end
        end
        iter = iter + 1; 
        if iter > max_iter
            break; 
        end
    end
     
end

