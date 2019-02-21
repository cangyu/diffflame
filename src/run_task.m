clear; close all; clc;

mkdir '../data' 'ChemTbl'

data = load('../data/task.txt');
n = length(data);

for k = 1: min(n, 500)
    %Clean environment
    delete('../data/iter*.txt');
	delete('../pic/*.png');
    
    %Run case
    mf = data(k, 1);
    mo = data(k, 2);
    try
        SolveFlame3(mf, mo, 0.05, 5001, '../data', 0);
    catch
        fprintf('Case%d failed!\n', k);
    end
end
