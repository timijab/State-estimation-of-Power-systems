%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN HELP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: [Vtrue Angletrue]= IEEE_true_value(num) %%
%%The function zconv returns the true system state  %%

function [Vtrue Angletrue]= IEEE_true_value(num)
IEEE_true_value_5=[1 1.010 0.00;
    2 0.995 -8.65;
    3 1.055 -10.75;
    4 1.090 -10.75;
    5 1.054 -11.81;];
IEEE_true_value_14= [1	1.06	0	2.323929	-0.165494	0	0	3	;
2	1.045	-4.98	0.4	0.435567	0.217	0.127	2	;
3	1.01	-12.72	0	0.250753	0.942	0.19	2	;
4	1.019	-10.32	0	0	0.478	-0.039	0	;
5	1.020	-8.779	0	0	0.076	0.016	0	;
6	1.070	-14.22	0	0.127298	0.112	0.075	2	;
7	1.062	-13.37	0	0	0	0	0	;
8	1.090	-13.37	0	0.17623	0	0	2	;
9	1.056	-14.95	0	0	0.295	0.166	0	;
10	1.051	-15.10	0	0	0.09	0.058	0	;
11	1.057	-14.79	0	0	0.035	0.018	0	;
12	1.055	-15.08	0	0	0.061	0.016	0	;
13	1.050	-15.16	0	0	0.135	0.058	0	;
14	1.036	-16.04	0	0	0.149	0.05	0	;];
%%%% IEEE 30 Bus System    %%%%
IEEE_true_value_30=[1	1.06	0	2.610866	-0.317445	0	0	3	;
2	1.043	-5.30	0.4	0.243661	0.217	0.127	2	;
3	1.021	-7.53	0	0	0.024	0.012	0	;
4	1.012	-9.28	0	0	0.076	0.016	0	;
5	1.010	-14.16	0	0.306673	0.942	0.19	2	;
6	1.010	-11.06	0	0	0	0	0	;
7	1.002	-12.86	0	0	0.228	0.109	0	;
8	1.010	-11.81	0	0.034272	0.3	0.3	2	;
9	1.051	-14.11	0	0	0	0	0	;
10	1.045	-15.70	0	0	0.058	0.02	0	;
11	1.082	-14.11	0	0.105345	0	0	2	;
12	1.057	-14.94	0	0	0.112	0.075	0	;
13	1.071	-14.94	0	-0.094343	0	0	2	;
14	1.042	-15.83	0	0	0.062	0.016	0	;
15	1.038	-15.92	0	0	0.082	0.025	0	;
16	1.044	-15.52	0	0	0.035	0.018	0	;
17	1.040	-15.86	0	0	0.09	0.058	0	;
18	1.028	-16.54	0	0	0.032	0.009	0	;
19	1.026	-16.71	0	0	0.095	0.034	0	;
20	1.030	-16.51	0	0	0.022	0.007	0	;
21	1.033	-16.14	0	0	0.175	0.112	0	;
22	1.033	-16.12	0	0	0	0	0	;
23	1.027	-16.31	0	0	0.032	0.016	0	;
24	1.022	-16.49	0	0	0.087	0.067	0	;
25	1.017	-16.07	0	0	0	0	0	;
26	1.000	-16.49	0	0	0.035	0.023	0	;
27	1.023	-15.55	0	0	0	0	0	;
28	1.007	-11.69	0	0	0	0	0	;
29	1.003	-16.78	0	0	0.024	0.009	0	;
30	0.992	-17.66	0	0	0.106	0.019	0	];
%-------------------------------------------------%
% True values the nigerian bus system%
Nig_true_56= [1	1.0279	-0.59284;
2	1.0068	-0.5105;
3	1.0269	-0.58602;
4	0.98272	-1.0896;
5	1.0018	-0.4906;
6	0.99541	-0.69624;
7	1.0408	-0.47789;
8	0.99257	-0.24255;
9	1.0213	-0.34787;
10	1.0236	-0.96643;
11	1.0081	-0.66555;
12	1.032	-0.56554;
13	1.0167	-0.98558;
14	1.0178	-0.98743;
15	0.99583	-0.20536;
16	1           0;
17	0.99189	-0.82507;
18	1       -0.46368;
19	1.0199	-0.79751;
20	1       -0.46475;
21	1.03	-0.78188;
22	1.0309	-0.80197;
23	0.99397	-0.21576;
24	1.0015	-0.49137;
25	1.0078	-0.4816;
26	1.03	-0.54647;
27	1       -0.55599;
28	1.0257	-0.58425;
29	0.99631	-0.23574;
30	1.0286   0.94768;
31	1.04	-0.91859;
32	1.0156	-0.95261;
33	1.0322	-1.008;
34	1.0267	-1.0241;
35	1       -0.63248;
36	1.03	-0.91404;
37	1.0112	-0.72648;
38	1       -0.68217;
39	1.03	-0.77721;
40	1.03	-0.23387;
41	0.9923	-0.2432;
42	0.99411	-0.83002;
43	0.98567	-94774;
44	1.0405	-0.81904;
45	1.0377	-0.82235;
46	1.03	-0.77189;
47	1       -0.47804;
48	1.0105	-0.50311;
49	1.03	-0.26716;
50	1.0269	-0.28879;
51	1.0479	-0.56047;
52	1.0189	-0.98936;
53	1.0188	-0.98996;
54	0.99458	-0.22349;
55	1.0279	-0.26074;
56	1.0443	-0.79833;];

%%%% IEEE 118 Bus System %%%%
% |Msnt | Type | Value | From | To | Rii | 
IEEE_true_value_118= [ 1    0.9550   -19.02;
   2    0.9714   -18.48;
   3    0.9677   -18.14;
   4    0.9980   -14.42;
   5    1.0020   -13.9722;
   6    0.9900   -16.70;
   7    0.9894   -17.1434;
   8    1.0150    -8.9520;
   9    1.043    -1.6984;
  10    1.0500     5.8821;
  11    0.9851   -16.99;
  12    0.9900   -17.5015;
  13    0.9683   -18.3600;
  14    0.984   -18.22;
  15    0.9700   -18.5113;
  16    0.9840   -17.8033;
  17    0.9952   -15.9953;
  18    0.9730   -18.2086;
  19    0.9635   -18.6932;
  20    0.9581   -17.8133;
  21    0.9587   -16.2257;
  22    0.9697   -13.6714;
  23    0.9997    -8.7511;
  24    0.9920    -8.8833;
  25    1.0500    -1.82;
  26    1.0150    -0.040;
  27    0.9680   -14.3887;
  28    0.9616   -16.1133;
  29    0.9633   -17.1058;
  30    0.9856   -10.9610;
  31    0.9670   -16.9890;
  32    0.9636   -14.9414;
  33    0.9716   -19.14;
  34    0.9859   -18.4942;
  35    0.9807   -18.9221;
  36    0.9801   -18.9176;
  37    0.9920   -18.0302;
  38    0.9620   -12.8938;
  39    0.9705   -21.4027;
  40    0.9700   -22.4764;
  41    0.9669   -22.9226;
  42    0.9850   -21.3270;
  43    0.9785   -18.5404;
  44    0.9851   -16.0547;
  45    0.9867   -14.2240;
  46    1.0050   -11.4182;
  47    1.0171    -9.1961;
  48    1.0207    -9.9758;
  49    1.0250    -8.9728;
  50    1.0011   -11.0117;
  51    0.9669   -13.6302;
  52    0.9569   -14.5834;
  53    0.9460   -15.5580;
  54    0.9550   -14.6459;
  55    0.9520   -14.9359;
  56    0.9540   -14.7499;
  57    0.9706   -13.5457;
  58    0.9591   -14.4021;
  59    0.9850   -10.5477;
  60    0.9932    -6.7670;
  61    0.9950    -5.8759;
  62    0.9980    -6.4918;
  63    0.9688    -7.1695;
  64    0.9838    -5.41;
  65    1.0050    -2.2794;
  66    1.0500    -2.4383;
  67    1.0197    -5.0778;
  68    1.0033    -2.4006;
  69    1.0350     0.0000;
  70    0.9840    -7.3808;
  71    0.9869    -7.7915;
  72    0.9800    -8.8888;
  73    0.9910    -8.0029;
  74    0.9580    -8.3295;
  75    0.9674    -7.0676;
  76    0.9430    -8.1969;
  77    1.0060    -3.2422;
  78    1.0034    -3.5461;
  79    1.0092    -3.2473;
  80    1.0400    -1.0025;
  81    0.9968    -1.8515;
  82    0.9888    -2.7219;
  83    0.9846    -1.5321;
  84    0.9798     1.0008;
  85    0.9850     2.56;
  86    0.9867     1.19;
  87    1.0150     1.45;
  88    0.9875     5.69;
  89    1.0050     9.7400;
  90    0.9850     3.3401;
  91    0.9800     3.3595;
  92    0.9924     3.8574;
  93    0.9870     0.8460;
  94    0.9906    -1.3082;
  95    0.9810    -2.2804;
  96    0.9927    -2.4481;
  97    1.0114    -2.0754;
  98    1.0235    -2.5528;
  99    1.0100    -2.9136;
 100    1.0170    -1.9167;
 101    0.9925    -0.3435;
 102    0.9911     2.36;
 103    1.0007    -5.5116;
 104    0.9710    -8.2563;
 105    0.9660    -9.3785;
 106    0.9618    -9.6264;
 107    0.9520   -12.4163;
 108    0.9668   -10.5694;
 109    0.968   -11.0186;
 110    0.9730   -11.8576;
 111    0.9800   -10.2126;
 112    0.9750   -14.9566;
 113    0.9930   -15.9973;
 114    0.9605   -15.2714;
 115    0.9600   -15.2790;
 116    1.0050    -2.8356;
 117    0.9739   -19.0424;
 118    0.9495    -8.0550;];
%-----------------------------------------------------%
%Return the data for the selected bus system
switch num
case 5
Vtrue=IEEE_true_value_5(:,2);
Angletrue=IEEE_true_value_5(:,3);
case 14
Vtrue=IEEE_true_value_14(:,2);
Angletrue=IEEE_true_value_14(:,3);
case 30
Vtrue=IEEE_true_value_30(:,2);
Angletrue=IEEE_true_value_30(:,3);
case 118
Vtrue=IEEE_true_value_118(:,2);
Angletrue=IEEE_true_value_118(:,3);
end
end