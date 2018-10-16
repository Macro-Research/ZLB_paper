function [regimeNew]= findRegime(regimeOld,p_LN,p_NL)

draw=rand;
if regimeOld==1;%if previous regime was 1 (L regime)
    if draw<(1-p_LN);%no transition with prob 1-p_LN
        regimeNew=1;
    else regimeNew=0;
    end
else %if previous regime was 0 (N regime)
    if draw<p_NL %transition to regime 1 with prob p_NL
        regimeNew=1;
    else regimeNew=0;
    end
end



end