clear; close all; clc;

mkdir '../data' 'ChemTbl'

data = load('../data/task.txt');
for k = 1:length(data) 
    %Clean environment
    delete('../data/*.txt');
	delete('../pic/*.png');
    
    %Run case
    mf = data(k, 1);
    mo = data(k, 2);
    try
        SolveFlame1(mf, mo, 0.05, 2501, '../data/ChemTbl');
    catch
        fprintf('Case%d failed!\n', k);
    end
end
