function[cm]=correction_factor(m)
if m>9
     if mod(m,2)==0  % the number is even ?
       cm=m/(m-0.9);
     else
         cm=1; % the number m is odd.
     end
    else
      switch m
        case 2;cm=0.743;
        case 3;cm=1.851;
        case 4;cm=0.954;
        case 5;cm=1.351;
        case 6;cm=0.993;
        case 7;cm=1.198;
        case 8;cm=1.005;
        case 9;cm=1.131;
      end
end
end