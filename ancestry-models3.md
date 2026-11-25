Ancestry Model for Classifying White River Spring-Run Natural Origin
Chinook Salmon
================
Tyler McCraney

#### White River NORs & White/Puyallup (Green River Origin) Baseline

    ## [1] 251

``` zsh
conda activate admixture_Env

cd /Users/tylermccraney/Google\ Drive/My\ Drive/Puyallup-White/admixture

plink --file gdatw3 --make-bed --out gdatw3

for K in 1 2 3 4; do admixture -s $RANDOM --cv=10 gdatw3.bed $K | tee gdatw3_log${K}.out; done

grep -h CV gdatw3_log* > gdatw3_log.txt
```

![](ancestry-models3_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

![](ancestry-models3_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Sort ancestry. The ‘Fall’ class is defined by 0-40% ancestry, the
‘Hybrid’ class by 40-60% ancestry, and the ‘Spring’ class by 60-100%
ancestry.

![](ancestry-models3_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Baseline reporting groups (rows) assigned to ‘Spring’, ‘Fall’, or
‘Hybrid’ ancestry classes (columns). The ‘Fall’ class is defined by
0-40% ancestry, the ‘Hybrid’ class by 40-60% ancestry, and the ‘Spring’
class by 60-100% ancestry.

    ##         
    ##          Spring Hybrid Fall
    ##   Spring    214     11    0
    ##   Fall        1     18  237

White River NOR assignments (rows) assigned to ‘Spring’, ‘Fall’, or
‘Hybrid’ ancestry classes (columns). Ancestry classes defined same as
above.

    ##          
    ##           Spring Hybrid Fall
    ##   Spring     712     29    0
    ##   Unknown     15     19    2
    ##   Fall        10     33   15
