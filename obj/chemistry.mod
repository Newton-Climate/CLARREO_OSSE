V24 chemistry
13 chemistry.F90 S579 0
02/18/2009  14:54:02
use abortutils private
use ghg_surfvals private
use constituents private
use physconst private
use ppgrid private
use pmgrid private
use shr_kind_mod private
enduse
D 125 21 9 3 262 261 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 189 26 3 189 189
D 128 21 9 3 262 261 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 189 26 3 189 189
D 131 21 9 3 262 261 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 189 26 3 189 189
D 134 21 9 3 262 261 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 189 26 3 189 189
D 137 21 9 3 262 23 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 36 26 3 36 36
D 140 21 9 3 262 23 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 36 26 3 36 36
D 143 21 9 3 262 23 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 36 26 3 36 36
D 146 21 9 3 262 23 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
 0 36 26 3 36 36
D 149 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 152 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 155 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 158 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 161 21 6 1 3 189 0 0 0 0 0
 0 189 3 3 189 189
D 164 21 6 1 3 189 0 0 0 0 0
 0 189 3 3 189 189
D 167 18 178
D 169 18 13
D 171 21 169 1 3 15 0 0 0 0 0
 0 15 3 3 15 15
D 174 21 167 1 3 15 0 0 0 0 0
 0 15 3 3 15 15
D 177 21 169 1 3 15 0 0 0 0 0
 0 15 3 3 15 15
S 579 24 0 0 0 8 1 0 4653 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 chemistry
S 581 23 0 0 0 8 618 579 4676 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 r8
S 584 23 0 0 0 8 638 579 4698 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 plat
S 585 23 0 0 0 8 637 579 4703 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 plev
S 586 23 0 0 0 8 641 579 4708 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 plevp
S 587 23 0 0 0 8 636 579 4714 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 plon
S 588 23 0 0 0 6 661 579 4719 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 masterproc
S 590 23 0 0 0 8 664 579 4737 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 pcols
S 591 23 0 0 0 8 665 579 4743 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 pver
S 593 23 0 0 0 6 790 579 4758 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwdry
S 594 23 0 0 0 6 794 579 4764 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwch4
S 595 23 0 0 0 6 793 579 4770 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwn2o
S 596 23 0 0 0 6 795 579 4776 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwf11
S 597 23 0 0 0 6 796 579 4782 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwf12
S 598 23 0 0 0 6 792 579 4788 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 mwh2o
S 599 23 0 0 0 8 777 579 4794 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 gravit
S 601 23 0 0 0 8 822 579 4814 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ppcnst
S 602 23 0 0 0 8 855 579 4821 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cnst_add
S 603 23 0 0 0 8 825 579 4830 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cnst_name
S 604 23 0 0 0 8 823 579 4840 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 advected
S 606 23 0 0 0 8 887 579 4862 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ch4vmr
S 607 23 0 0 0 6 886 579 4869 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 n2ovmr
S 608 23 0 0 0 8 888 579 4876 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 f11vmr
S 609 23 0 0 0 8 889 579 4883 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 f12vmr
S 611 23 0 0 0 8 807 579 4901 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 endrun
S 613 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 614 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 616 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 80 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 618 16 2 shr_kind_mod shr_kind_r8
S 626 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 630 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 636 16 1 pmgrid plon
R 637 16 2 pmgrid plev
R 638 16 3 pmgrid plat
R 641 16 6 pmgrid plevp
R 661 16 26 pmgrid masterproc
R 664 16 1 ppgrid pcols
R 665 16 2 ppgrid pver
R 777 16 20 physconst gravit
R 790 16 33 physconst mwdry
R 792 16 35 physconst mwh2o
R 793 16 36 physconst mwn2o
R 794 16 37 physconst mwch4
R 795 16 38 physconst mwf11
R 796 16 39 physconst mwf12
R 807 14 5 abortutils endrun
S 808 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 822 16 13 constituents ppcnst
R 823 16 14 constituents advected
R 825 7 16 constituents cnst_name
R 855 14 46 constituents cnst_add
S 876 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 886 6 9 ghg_surfvals n2ovmr
R 887 6 10 ghg_surfvals ch4vmr
R 888 6 11 ghg_surfvals f11vmr
R 889 6 12 ghg_surfvals f12vmr
S 971 27 0 0 0 8 1034 579 7072 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 chem_register
S 972 6 4 0 0 16 1 579 7086 80000c 0 0 0 0 0 0 0 0 0 1029 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 trace_gas
S 973 16 0 0 0 6 1 579 7096 14 400000 0 0 0 0 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ptrlon
S 974 16 0 0 0 6 1 579 7103 14 400000 0 0 0 0 36 246 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ptrlat
S 975 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 36 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 976 16 0 0 0 6 1 579 7110 14 400000 0 0 0 0 56 248 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ptrlev
S 977 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 56 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 978 16 1 0 0 6 1 579 7117 14 400000 0 0 0 0 12 189 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ptrtim
S 979 16 0 0 0 9 1 579 7124 14 400000 0 0 0 0 980 251 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 rmwn2o
S 980 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1073237482 1414554581 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 981 16 0 0 0 9 1 579 7131 14 400000 0 0 0 0 982 253 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 rmwch4
S 982 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1071754503 -2094845611 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 983 16 0 0 0 9 1 579 7138 14 400000 0 0 0 0 984 255 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 rmwf11
S 984 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1074972631 -78289814 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 985 16 0 0 0 9 1 579 7145 14 400000 0 0 0 0 986 257 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 rmwf12
S 986 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1074827831 183565888 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 987 16 0 0 0 9 1 579 7152 14 400000 0 0 0 0 988 259 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 rh2och4
S 988 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1072825368 -1821066134 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 989 7 4 0 4 125 992 579 7160 800014 100 0 0 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tch4i
S 990 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 480 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 991 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 42 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 992 7 4 0 4 128 993 579 7166 800014 100 0 3840 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tn2oi
S 993 7 4 0 4 131 994 579 7172 800014 100 0 7680 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc11i
S 994 7 4 0 4 134 995 579 7180 800014 100 0 11520 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc12i
S 995 7 4 0 4 137 996 579 7188 800014 100 0 15360 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tch4m
S 996 7 4 0 4 140 997 579 7194 800014 100 0 16000 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tn2om
S 997 7 4 0 4 143 998 579 7200 800014 100 0 16640 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc11m
S 998 7 4 0 4 146 999 579 7208 800014 100 0 17280 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc12m
S 999 7 4 0 4 149 1000 579 7216 800014 100 0 17920 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tch4
S 1000 7 4 0 4 152 1001 579 7221 800014 100 0 18240 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tn2o
S 1001 7 4 0 4 155 1002 579 7226 800014 100 0 18560 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc11
S 1002 7 4 0 4 158 1003 579 7233 800014 100 0 18880 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 tcfc12
S 1003 6 4 0 0 9 1004 579 7240 14 0 0 19200 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cdaytrm
S 1004 6 4 0 0 9 1 579 7248 14 0 0 19208 0 0 0 0 0 0 1030 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cdaytrp
S 1005 6 4 0 0 6 1006 579 7256 14 0 0 0 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 np
S 1006 6 4 0 0 6 1007 579 7259 14 0 0 4 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 nm
S 1007 6 4 0 0 6 1008 579 7262 14 0 0 8 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 np1
S 1008 7 4 0 4 161 1009 579 7266 800014 100 0 16 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 date_tr
S 1009 7 4 0 4 164 1028 579 7274 800014 100 0 64 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 sec_tr
S 1010 16 0 0 0 9 1 579 7281 14 400000 0 0 0 0 1012 264 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cpch4
S 1012 3 0 0 0 9 0 1 0 0 0 0 0 0 0 1082445824 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 1013 16 0 0 0 9 1 579 7292 14 400000 0 0 0 0 1012 264 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cpn2o
S 1014 16 0 0 0 9 1 579 7303 14 400000 0 0 0 0 1012 264 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cpf11
S 1015 16 0 0 0 9 1 579 7314 14 400000 0 0 0 0 1012 264 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cpf12
S 1016 16 1 0 0 6 1 579 7325 14 400000 0 0 0 0 4 15 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ncnst
S 1021 16 0 0 0 171 1 579 7355 14 400000 0 0 0 0 1026 1 0 0 0 0 0 0 0 0 0 0 0 282 0 579 0 0 0 0 cnst_names
S 1022 3 0 0 0 169 0 1 0 0 0 0 0 0 0 7366 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 8 4e 32 4f 20 20 20 20 20
S 1023 3 0 0 0 169 0 1 0 0 0 0 0 0 0 7375 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 8 43 48 34 20 20 20 20 20
S 1024 3 0 0 0 169 0 1 0 0 0 0 0 0 0 7384 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 8 43 46 43 31 31 20 20 20
S 1025 3 0 0 0 169 0 1 0 0 0 0 0 0 0 7393 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 8 43 46 43 31 32 20 20 20
S 1026 7 4 0 4 171 1 579 7402 4080005c 400100 0 0 0 0 0 0 0 0 1032 0 0 0 0 0 0 0 0 282 0 579 0 1021 0 0 cnst_names$ac
S 1027 7 4 0 4 177 1 579 7416 800014 100 0 0 0 0 0 0 0 0 1033 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 srcnam
S 1028 6 4 0 0 6 1 579 7423 14 0 0 112 0 0 0 0 0 0 1031 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 ixghg
S 1029 11 0 0 0 8 954 579 7429 40800000 801000 0 4 0 0 972 972 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chemistry$8
S 1030 11 0 0 1 8 1029 579 7441 40800010 801000 0 19216 0 0 989 1004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chemistry$6
S 1031 11 0 0 1 8 1030 579 7453 40800010 801000 0 116 0 0 1005 1028 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chemistry$4
S 1032 11 0 0 1 8 1031 579 7465 40800010 801000 0 32 0 0 1026 1026 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chemistry$13
S 1033 11 0 0 1 8 1032 579 7478 40800010 801000 0 32 0 0 1027 1027 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chemistry$5
S 1034 23 5 0 0 0 1035 579 7072 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 chem_register
S 1035 14 5 0 0 0 1 1034 7072 0 400000 0 0 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 90 0 579 0 0 0 0 chem_register
F 1035 0
A 13 2 0 0 0 6 613 0 0 0 13 0 0 0 0 0 0 0 0 0
A 15 2 0 0 0 6 614 0 0 0 15 0 0 0 0 0 0 0 0 0
A 23 2 0 0 0 6 616 0 0 0 23 0 0 0 0 0 0 0 0 0
A 26 2 0 0 0 6 626 0 0 0 26 0 0 0 0 0 0 0 0 0
A 36 2 0 0 0 6 630 0 0 0 36 0 0 0 0 0 0 0 0 0
A 178 2 0 0 0 6 808 0 0 0 178 0 0 0 0 0 0 0 0 0
A 189 2 0 0 0 6 876 0 0 0 189 0 0 0 0 0 0 0 0 0
A 246 2 0 0 0 6 975 0 0 0 246 0 0 0 0 0 0 0 0 0
A 248 2 0 0 204 6 977 0 0 0 248 0 0 0 0 0 0 0 0 0
A 251 2 0 0 0 9 980 0 0 0 251 0 0 0 0 0 0 0 0 0
A 253 2 0 0 0 9 982 0 0 0 253 0 0 0 0 0 0 0 0 0
A 255 2 0 0 0 9 984 0 0 0 255 0 0 0 0 0 0 0 0 0
A 257 2 0 0 0 9 986 0 0 0 257 0 0 0 0 0 0 0 0 0
A 259 2 0 0 0 9 988 0 0 0 259 0 0 0 0 0 0 0 0 0
A 261 2 0 0 0 6 990 0 0 0 261 0 0 0 0 0 0 0 0 0
A 262 2 0 0 0 6 991 0 0 0 262 0 0 0 0 0 0 0 0 0
A 264 2 0 0 17 9 1012 0 0 0 264 0 0 0 0 0 0 0 0 0
A 274 2 0 0 200 169 1022 0 0 0 274 0 0 0 0 0 0 0 0 0
A 275 15 0 0 0 167 1021 274 277 0 0 0 0 0 0 0 0 0 0 0
A 276 2 0 0 176 169 1023 0 0 0 276 0 0 0 0 0 0 0 0 0
A 277 15 0 0 0 167 1021 276 279 0 0 0 0 0 0 0 0 0 0 0
A 278 2 0 0 102 169 1024 0 0 0 278 0 0 0 0 0 0 0 0 0
A 279 15 0 0 0 167 1021 278 281 0 0 0 0 0 0 0 0 0 0 0
A 280 2 0 0 104 169 1025 0 0 0 280 0 0 0 0 0 0 0 0 0
A 281 15 0 0 0 167 1021 280 0 0 0 0 0 0 0 0 0 0 0 0
A 282 15 0 0 0 174 1021 275 0 0 0 0 0 0 0 0 0 0 0 0
A 283 1 0 3 25 171 1026 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 82 1 1
V 283 171 7 0
R 0 174 0 0
A 0 169 0 0 1 274 1
A 0 169 0 0 1 276 1
A 0 169 0 0 1 278 1
A 0 169 0 0 1 280 0
Z
