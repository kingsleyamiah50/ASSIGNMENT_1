function [outputArg1] = myRand(nElectrons)


    
    Vtherm = 1.87e05;
    
    myNum(1:nElectrons) = randn(nElectrons,1) * Vtherm;
    
    myResult(1:nElectrons) = myNum(1: nElectrons) * -1;
    for i = 1:nElectrons
        if(myResult(i) == abs(myNum(i)))
            outputArg1(i) = -1;
        else
            outputArg1(i) = 1;
        end
    end
end

