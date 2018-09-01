function [regimeNew]= findRegime(regimeOld,p_11,p_22)

draw=rand;
if regimeOld==1;%if previous regime was 1 (L regime)
    if draw<(p_11);%no transition with prob p_11
        regimeNew=1;
    else regimeNew=0;
    end
else %if previous regime was 0 (N regime)
    if draw<p_22 %no transition with prob p_22
        regimeNew=0;
    else regimeNew=1;
    end
end



end