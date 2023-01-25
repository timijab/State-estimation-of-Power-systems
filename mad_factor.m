function[bm]=mad_factor(m)
if m>9
    bm=m/(m-0.8);
else
    switch m
        case 1;bm=1.000;
        case 2;bm=1.196;
        case 3;bm=1.495;
        case 4;bm=1.363;
        case 5;bm=1.206;
        case 6;bm=1.200;
        case 7;bm=1.140;
        case 8;bm=1.129;
        case 9;bm=1.107;
    end
end