V24 shr_sys_mod
15 shr_sys_mod.F90 S579 0
02/18/2009  12:19:42
use shr_kind_mod private
enduse
S 579 24 0 0 0 8 1 0 4653 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 shr_sys_mod
S 581 23 0 0 0 8 596 579 4678 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_kind_in
S 582 23 0 0 0 8 591 579 4690 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_kind_r8
S 583 23 0 0 0 8 594 579 4702 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_kind_i8
R 591 16 2 shr_kind_mod shr_kind_r8
R 594 16 5 shr_kind_mod shr_kind_i8
R 596 16 7 shr_kind_mod shr_kind_in
S 979 27 0 0 0 8 986 579 6877 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_system
S 980 27 0 0 0 8 990 579 6892 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_chdir
S 981 27 0 0 0 8 994 579 6906 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_getenv
S 982 27 0 0 0 8 999 579 6921 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_flush
S 983 27 0 0 0 8 1002 579 6935 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_abort
S 984 27 0 0 0 8 1006 579 6949 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_irtc
S 985 27 0 0 0 8 1010 579 6962 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 shr_sys_sleep
S 986 23 5 0 0 0 989 579 6877 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_system
S 987 1 3 1 0 28 1 986 6976 14 43000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 str
S 988 1 3 2 0 6 1 986 6980 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rcode
S 989 14 5 0 0 0 1 986 6877 0 400000 0 0 0 174 2 0 0 0 0 0 0 0 0 0 0 0 0 26 0 579 0 0 0 0 shr_sys_system
F 989 2 987 988
S 990 23 5 0 0 0 993 579 6892 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_chdir
S 991 1 3 1 0 28 1 990 6986 14 43000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 path
S 992 1 3 2 0 6 1 990 6991 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rcode
S 993 14 5 0 0 0 1 990 6892 0 400000 0 0 0 177 2 0 0 0 0 0 0 0 0 0 0 0 0 73 0 579 0 0 0 0 shr_sys_chdir
F 993 2 991 992
S 994 23 5 0 0 0 998 579 6906 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_getenv
S 995 1 3 1 0 28 1 994 6997 14 43000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 name
S 996 1 3 2 0 28 1 994 7002 14 43000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 val
S 997 1 3 2 0 6 1 994 7006 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rcode
S 998 14 5 0 0 0 1 994 6906 0 400000 0 0 0 180 3 0 0 0 0 0 0 0 0 0 0 0 0 114 0 579 0 0 0 0 shr_sys_getenv
F 998 3 995 996 997
S 999 23 5 0 0 0 1001 579 6921 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_flush
S 1000 1 3 0 0 6 1 999 3876 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 unit
S 1001 14 5 0 0 0 1 999 6921 0 400000 0 0 0 184 1 0 0 0 0 0 0 0 0 0 0 0 0 155 0 579 0 0 0 0 shr_sys_flush
F 1001 1 1000
S 1002 23 5 0 0 0 1005 579 6935 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_abort
S 1003 1 3 0 0 28 1 1002 7012 80000014 43000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 string
S 1004 1 3 0 0 6 1 1002 7019 80000014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rcode
S 1005 14 5 0 0 0 1 1002 6935 0 400000 0 0 0 186 2 0 0 0 0 0 0 0 0 0 0 0 0 181 0 579 0 0 0 0 shr_sys_abort
F 1005 2 1003 1004
S 1006 23 5 0 0 7 1008 579 6949 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_irtc
S 1007 1 3 0 0 7 1 1006 7025 80000014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rate
S 1008 14 5 0 0 7 1 1006 6949 4 400000 0 0 0 189 1 0 0 1009 0 0 0 0 0 0 0 0 0 223 0 579 0 0 0 0 shr_sys_irtc
F 1008 1 1007
S 1009 1 3 0 0 7 1 1006 6949 14 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_irtc
S 1010 23 5 0 0 0 1012 579 6962 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shr_sys_sleep
S 1011 1 3 1 0 9 1 1010 7030 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sec
S 1012 14 5 0 0 0 1 1010 6962 0 400000 0 0 0 191 1 0 0 0 0 0 0 0 0 0 0 0 0 256 0 579 0 0 0 0 shr_sys_sleep
F 1012 1 1011
Z
Z
