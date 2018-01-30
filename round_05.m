function tn = round_05 (tt)

    if tt>0

       if (tt-fix(tt))>0.5

           if (tt-fix(tt)-0.5)/2>0.125

               tn=fix(tt)+1;

           else

               tn=fix(tt)+.5;

           end

       else

           if (tt-fix(tt))/2>0.125

               tn=fix(tt)+.5;

           else

               tn=fix(tt);

           end

       end

    elseif tt == 0
        tn = tt;
    end
end