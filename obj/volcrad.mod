V24 volcrad
11 volcrad.F90 S579 0
02/18/2009  14:54:00
use shr_kind_mod private
enduse
D 33 21 9 1 3 43 0 0 0 0 0
 0 43 3 3 43 43
D 39 21 9 2 46 29 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
D 42 21 9 4 64 66 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
 0 29 29 3 29 29
 0 43 65 3 43 43
D 45 21 9 2 46 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 48 21 9 2 46 29 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
S 579 24 0 0 0 8 1 0 4653 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 volcrad
S 581 23 0 0 0 8 591 579 4674 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 r8
S 587 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 591 16 2 shr_kind_mod shr_kind_r8
S 599 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 600 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 41 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 620 27 0 0 0 8 644 579 4895 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 aer_trn
S 621 27 0 0 0 8 651 579 4903 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 aer_pth
S 622 16 1 0 0 6 1 579 4911 4 400000 0 0 0 0 7 43 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 bnd_nbr_lw
S 623 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 624 16 0 0 0 6 1 579 4922 4 400000 0 0 0 0 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_h2o_nonwnd
S 625 16 0 0 0 6 1 579 4940 4 400000 0 0 0 0 2 46 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_h2o_window
S 626 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 627 16 0 0 0 6 1 579 4958 4 400000 0 0 0 0 3 48 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_0500_0650
S 628 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 629 16 0 0 0 6 1 579 4975 4 400000 0 0 0 0 4 15 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_0650_0800
S 630 16 0 0 0 6 1 579 4992 4 400000 0 0 0 0 5 51 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_0800_1000
S 631 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 632 16 0 0 0 6 1 579 5009 4 400000 0 0 0 0 6 53 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_1000_1200
S 633 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 634 16 0 0 0 6 1 579 5026 4 400000 0 0 0 0 7 43 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 idx_lw_1200_2000
S 635 7 4 0 4 33 1 579 5043 80000c 100 0 0 0 0 0 0 0 0 643 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 abs_cff_mss_aer
S 643 11 0 0 1 8 611 579 5150 40800000 801000 0 56 0 0 635 635 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 volcrad$10
S 644 23 5 0 0 0 647 579 4895 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_trn
S 645 7 3 1 0 39 1 644 5161 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_mpp
S 646 7 3 2 0 42 1 644 5169 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_trn_ttl
S 647 14 5 0 0 0 1 644 4895 0 400000 0 0 0 2 2 0 0 0 0 0 0 0 0 0 0 0 0 82 0 579 0 0 0 0 aer_trn
F 647 2 645 646
S 648 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 1681 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 649 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 11767 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 650 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 1724 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 651 23 5 0 0 0 655 579 4903 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_pth
S 652 7 3 1 0 45 1 651 5181 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_mass
S 653 7 3 2 0 48 1 651 5190 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 aer_mpp
S 654 1 3 1 0 6 1 651 5198 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 655 14 5 0 0 0 1 651 4903 0 400000 0 0 0 5 3 0 0 0 0 0 0 0 0 0 0 0 0 162 0 579 0 0 0 0 aer_pth
F 655 3 652 653 654
A 15 2 0 0 0 6 587 0 0 0 15 0 0 0 0 0 0 0 0 0
A 26 2 0 0 0 6 599 0 0 0 26 0 0 0 0 0 0 0 0 0
A 29 2 0 0 0 6 600 0 0 0 29 0 0 0 0 0 0 0 0 0
A 43 2 0 0 0 6 623 0 0 0 43 0 0 0 0 0 0 0 0 0
A 46 2 0 0 0 6 626 0 0 0 46 0 0 0 0 0 0 0 0 0
A 48 2 0 0 0 6 628 0 0 0 48 0 0 0 0 0 0 0 0 0
A 51 2 0 0 0 6 631 0 0 0 51 0 0 0 0 0 0 0 0 0
A 53 2 0 0 0 6 633 0 0 0 53 0 0 0 0 0 0 0 0 0
A 64 2 0 0 0 6 650 0 0 0 64 0 0 0 0 0 0 0 0 0
A 65 2 0 0 0 6 648 0 0 0 65 0 0 0 0 0 0 0 0 0
A 66 2 0 0 0 6 649 0 0 0 66 0 0 0 0 0 0 0 0 0
Z
Z