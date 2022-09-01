function scaleSTAtemp=scaleSTA(STAtr,STAtempsplot,scale0)
%   scaleSItemp=scaleSI(SILBtr,SILBtempsplot,localmaxSI/maxSILBtemp);
scaleSTAtemp=fminunc(@nested,scale0);
    function y=nested(a)
        y=sum((a*STAtempsplot-STAtr').^2); %(160,1) vs (1,160).  And later not. And later yes.
        %y=sum(abs(a*STAtempsplot-STAtr')); %(160,1) vs (1,160).
    end
end