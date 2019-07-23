C = 1e3;
S = 10;
label = 5;
fp0 = zeros(C, S);
varfp0 = zeros(C, S);
for i = 1:C
    load([num2str(label) '/C' num2str(i) '/newborn_selected'])
    for j = 1:S
        fp0(i,j) = sum(newborn_selected(j).M_L.*newborn_selected(j).fp)/...
            sum(newborn_selected(j).M_L);
        varfp0(i,j) = sum(newborn_selected(j).M_L.*(newborn_selected(j).fp-fp0(i,j)).^2)/...
            sum(newborn_selected(j).M_L);
    end
end
%%
save([num2str(label) '/fp_data'], 'fp0', 'varfp0')
        
    