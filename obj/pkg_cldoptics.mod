V24 pkg_cldoptics
17 pkg_cldoptics.F90 S579 0
02/18/2009  14:54:00
use ppgrid private
use shr_kind_mod private
enduse
D 33 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 36 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 39 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 42 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 45 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 48 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 51 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 54 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 57 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 60 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 63 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 66 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 69 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 72 21 9 2 36 29 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
D 75 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 78 21 9 2 36 29 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
D 81 21 6 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 84 21 9 2 36 29 0 0 0 0 0
 0 3 3 3 3 3
 0 29 3 3 29 29
D 87 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 90 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 93 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 96 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 99 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 102 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 105 21 9 1 3 3 0 0 0 0 0
 0 3 3 3 3 3
D 108 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 111 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 114 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
D 117 21 9 2 36 26 0 0 0 0 0
 0 3 3 3 3 3
 0 26 3 3 26 26
S 579 24 0 0 0 8 1 0 4653 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 pkg_cldoptics
S 581 23 0 0 0 8 593 579 4680 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 r8
S 584 23 0 0 0 8 605 579 4702 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 pcols
S 585 23 0 0 0 8 606 579 4708 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 pver
S 586 23 0 0 0 8 608 579 4713 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 pverp
R 593 16 2 shr_kind_mod shr_kind_r8
S 601 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 602 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 41 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 605 16 1 ppgrid pcols
R 606 16 2 ppgrid pver
R 608 16 4 ppgrid pverp
S 614 6 4 0 0 8 615 579 4869 0 1000 0 0 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cldefr
S 615 6 4 0 0 8 616 579 4876 0 1000 0 4 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cldems
S 616 6 4 0 0 8 617 579 4883 0 1000 0 8 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cldovrlap
S 617 6 4 0 0 8 618 579 4893 0 1000 0 12 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 cldclw
S 618 6 4 0 0 8 619 579 4900 0 1000 0 16 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 reitab
S 619 6 4 0 0 8 1 579 4907 0 1000 0 20 0 0 0 0 0 0 620 0 0 0 0 0 0 0 0 0 0 579 0 0 0 0 reltab
S 620 11 0 0 0 8 613 579 4914 40800000 801000 0 24 0 0 614 619 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pkg_cldoptics$0
S 621 23 5 0 0 0 633 579 4869 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cldefr
S 622 1 3 1 0 6 1 621 4930 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 lchnk
S 623 1 3 1 0 6 1 621 4936 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 624 7 3 1 0 33 1 621 4941 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 landfrac
S 625 7 3 1 0 39 1 621 4950 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 626 7 3 2 0 54 1 621 4952 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rel
S 627 7 3 2 0 57 1 621 4956 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rei
S 628 7 3 1 0 42 1 621 4960 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ps
S 629 7 3 1 0 45 1 621 4963 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pmid
S 630 7 3 1 0 48 1 621 4968 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 landm
S 631 7 3 1 0 36 1 621 4974 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 icefrac
S 632 7 3 1 0 51 1 621 4982 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 snowh
S 633 14 5 0 0 0 1 621 4869 0 400000 0 0 0 2 11 0 0 0 0 0 0 0 0 0 0 0 0 21 0 579 0 0 0 0 cldefr
F 633 11 622 623 624 625 626 627 628 629 630 631 632
S 634 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 635 23 5 0 0 0 642 579 4876 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cldems
S 636 1 3 1 0 6 1 635 4988 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 lchnk
S 637 1 3 1 0 6 1 635 4994 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 638 7 3 1 0 60 1 635 4999 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 clwp
S 639 7 3 1 0 66 1 635 5004 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fice
S 640 7 3 1 0 63 1 635 5009 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rei
S 641 7 3 2 0 69 1 635 5013 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 emis
S 642 14 5 0 0 0 1 635 4876 0 400000 0 0 0 14 6 0 0 0 0 0 0 0 0 0 0 0 0 70 0 579 0 0 0 0 cldems
F 642 6 636 637 638 639 640 641
S 643 23 5 0 0 0 650 579 4883 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cldovrlap
S 644 1 3 1 0 6 1 643 5018 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 lchnk
S 645 1 3 1 0 6 1 643 5024 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 646 7 3 1 0 72 1 643 5029 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pint
S 647 7 3 1 0 75 1 643 5034 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cld
S 648 7 3 0 0 81 1 643 5038 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nmxrgn
S 649 7 3 2 0 78 1 643 5045 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pmxrgn
S 650 14 5 0 0 0 1 643 4883 0 400000 0 0 0 21 6 0 0 0 0 0 0 0 0 0 0 0 0 125 0 579 0 0 0 0 cldovrlap
F 650 6 644 645 646 647 648 649
S 651 23 5 0 0 0 658 579 4893 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cldclw
S 652 1 3 1 0 6 1 651 5052 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 lchnk
S 653 1 3 1 0 6 1 651 5058 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 654 7 3 1 0 84 1 651 5063 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 zi
S 655 7 3 0 0 90 1 651 5066 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 clwp
S 656 7 3 1 0 87 1 651 5071 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tpw
S 657 7 3 0 0 93 1 651 5075 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 hl
S 658 14 5 0 0 0 1 651 4893 0 400000 0 0 0 28 6 0 0 0 0 0 0 0 0 0 0 0 0 202 0 579 0 0 0 0 cldclw
F 658 6 652 653 654 655 656 657
S 659 23 5 0 0 0 667 579 4907 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 reltab
S 660 1 3 1 0 6 1 659 5078 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 661 7 3 1 0 108 1 659 5083 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 662 7 3 1 0 96 1 659 5085 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 landfrac
S 663 7 3 1 0 105 1 659 5094 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 landm
S 664 7 3 1 0 99 1 659 5100 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 icefrac
S 665 7 3 2 0 111 1 659 5108 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rel
S 666 7 3 1 0 102 1 659 5112 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 snowh
S 667 14 5 0 0 0 1 659 4907 0 400000 0 0 0 35 7 0 0 0 0 0 0 0 0 0 0 0 0 271 0 579 0 0 0 0 reltab
F 667 7 660 661 662 663 664 665 666
S 668 23 5 0 0 0 672 579 4900 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 reitab
S 669 1 3 1 0 6 1 668 5118 14 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ncol
S 670 7 3 1 0 117 1 668 5123 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 671 7 3 2 0 114 1 668 5125 800014 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 re
S 672 14 5 0 0 0 1 668 4900 0 400000 0 0 0 43 3 0 0 0 0 0 0 0 0 0 0 0 0 331 0 579 0 0 0 0 reitab
F 672 3 669 670 671
A 26 2 0 0 0 6 601 0 0 0 26 0 0 0 0 0 0 0 0 0
A 29 2 0 0 0 6 602 0 0 0 29 0 0 0 0 0 0 0 0 0
A 36 2 0 0 0 6 634 0 0 0 36 0 0 0 0 0 0 0 0 0
Z
Z