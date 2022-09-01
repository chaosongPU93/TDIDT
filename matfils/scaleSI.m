function scaleSItemp=scaleSI(SILBtr,SILBtempsplot,scale0)
%   scaleSItemp=scaleSI(SILBtr,SILBtempsplot,localmaxSI/maxSILBtemp);
scaleSItemp=fminunc(@nested,scale0);
    function y=nested(a)
        y=sum((a*SILBtempsplot-SILBtr').^2); %(160,1) vs (1,160)
    end
end