*-------------------------------------------------------------;
* Generate Linkage Term Project Problem for 621 Comp Stat Gen ;
* Author: M.A. Province 2005 (enhanced 2010)                  ;
*-------------------------------------------------------------;

options mprint mtrace;
title 'MLE Parametric Linkage Data';

                                                                                                                       
****************************************************************;
* Specify parameters                                           *;
****************************************************************;
*-----------------------------;
* Starting seed for simulation;
*-----------------------------;
%let seed=12345;
                                                                                                                       
*-----------------------------------------------;
* # Families (assume nuclear: F, M, C1, ..., Cm);
* # Kids/family                                 ;
*-----------------------------------------------;
%let nfams=100;
%let mkids=0;
                                                                                                                       
*--------------------------------------;
* Disease Gene (Latent)                ;
* Allele Frequency (assume BI-alllelic);
*--------------------------------------;
%let dq=0.4;
                                                                                                                       
*----------------------;
* Penetrance           ;
* Three genotypic means;
*----------------------;
%let mudaa=-0.6;
%let mudab=0.0;
%let mudbb=0.8;
                                                                                                                       
*----------------;
* Common variance;
*----------------;
%let de=0.1;
                                                                                                                       
*---------------------------------------;
* Linkage Marker (Observed)             ;
*   Number of Alleles for Linkage Marker;
*---------------------------------------;
%let lalleles=6;
                                                                                                                       
*------------------------------------------------------------------------------------------------------;
*Recombination fraction between Disease gene and Linkage Marker gene: 0=completely linked, 1/2=unlinked;
*------------------------------------------------------------------------------------------------------;
%let theta=0.1;
                                                                                                                       
****************************************************************;
* Done specifying parameters                                   *;
****************************************************************;
                                                                                                                       
                                                                                                                       
*----------------------------------------------------------------------------;
* Calculate Marker Allele Frequencies (assume all equal) from input # alleles;
*----------------------------------------------------------------------------;
data _null_;
  lafreq=(1/&lalleles);
  length lqs $200;
  lqs=trim(left(lafreq))||repeat(","||trim(left(lafreq)), &lalleles-2);
  call symput("lqs",trim(left(lqs)));
run;
                                                                                                                       
%put lqs=&lqs;
                                                                                                                       
*-----------------------------------------------------------------------------------------------------------;
* Macro LINK:                                                                                               ;
* Purpose:  Generate Linkage Dataset (Nuclear Families, Quantitative Pheno, Disease Gene and Marker)        ;
*-----------------------------------------------------------------------------------------------------------;
                                                                                                                       
%macro link;
                                                                                                                       
    title2 "theta=&theta mudaa=&mudaa mudab=&mudab mudbb=&mudbb dq=&dq de=&de lalleles=&lalleles nfams=&nfams mkids=&mkids seed=&seed";
                                                                                                                       
    *--------------------------------------------------------;
    * Generate &nfams Nuclear Families: F, M, 1, ..., &mkids.;
    * One record/pedid                                       ;
    *--------------------------------------------------------;
    data families;
      seed=&seed;
      dq=&dq;  * disease-risk allele freq of disease gene;
      dp=1-dq; * normal-risk allele freq of disease gene;
      de=&de;
                                                                                                                         
      *-------------------------------------------;
      * Disease Genotype, phenotype and haplotypes;
      *-------------------------------------------;
                                                                                                                         
      array dgene{*} dgf dgm dg1-dg&mkids;
      array phenos{*} pf  pm  p1 -p&mkids;
                                                                                                                         
      array mud{3} _temporary_ (&mudaa &mudab &mudbb);
                                                                                                                         
      *------------------------------------------------------------------;
      * Two alleles (blue from father, pink from mother) for Disease Gene;
      *------------------------------------------------------------------;
      array blue_dallele{*} blue_daf blue_dam blue_da1-blue_da&mkids;
      array pink_dallele{*} pink_daf pink_dam pink_da1-pink_da&mkids;
                                                                                                                         
      *--------------------------------------------------------------------;
      * Two alleles (blue from father, pink from mother) for Linkage Marker;
      *--------------------------------------------------------------------;
      array blue_lallele{*} blue_laf blue_lam blue_la1-blue_la&mkids;
      array pink_lallele{*} pink_laf pink_lam pink_la1-pink_la&mkids;
                                                                                                                         
      *--------------------;
      * Generate each pedid;
      *--------------------;
      do j=1 to &nfams;
                                                                                                                         
        pedid=j;
                                                                                                                         
        do i=1 to &mkids+2;
                                                                                                                         
          if i<=2 then do;
            *-------------------------------------------------------------------;
            * First generate Parents (founders) Genotypes--assume RANDOM mating ;
            *-------------------------------------------------------------------;
                                                                                                                         
            *---------------------------------------------------------;
            * Generate Disease Genotypes by drawing alleles (randomly);
            *---------------------------------------------------------;
            call rantbl(seed, dp, dq, blue_dallele{i});
            call rantbl(seed, dp, dq, pink_dallele{i});
                                                                                                                         
            *-------------------------------------------;
            * Generate Linkage Marker Alleles (randomly);
            *-------------------------------------------;
            call rantbl(seed, &lqs, blue_lallele{i});
            call rantbl(seed, &lqs, pink_lallele{i});

          end;
          else do;
            *----------------------------------------------------------;
            * Now generate Kids Genotypes from parents (dropping genes);
            *----------------------------------------------------------;
                                                                                                                         
            *-----------------------------------------------------------------------------------------------------------------;
            * Randomly choose whether father gives his Blue or Pink allele to each kid (this will be kids Blue allele)        ;
            * But you must do this for both the disease and linkage markers simultaneously, because they are linked by &theta ;
            *-----------------------------------------------------------------------------------------------------------------;
                                                                                                                        
            call ranuni(seed, xr1);
            blue_recomb=(xr1 < &theta);
            call ranuni(seed, x1);
                                                                                                                         
            if (x1<0.5) then do;
              blue_dallele{i} = blue_dallele{1};
              if blue_recomb then blue_lallele{i} = pink_lallele{1};
              else                blue_lallele{i} = blue_lallele{1};
            end;
            else do;
              blue_dallele{i} = pink_dallele{1};
              if blue_recomb then blue_lallele{i} = blue_lallele{1};
              else                blue_lallele{i} = pink_lallele{1};
            end;
                                                                                                                         
            *----------------------------------------;
            * Same for mother (pick kids Pink allele);
            *----------------------------------------;
                                                                                                                         
            call ranuni(seed, xr2);
            pink_recomb=(xr2 < &theta);
            call ranuni(seed, x2);
                                                                                                                         
            if (x2<0.5) then do;
              pink_dallele{i} = blue_dallele{2};
              if pink_recomb then pink_lallele{i} = pink_lallele{2};
              else                pink_lallele{i} = blue_lallele{2};
            end;
            else do;
              pink_dallele{i} = pink_dallele{2};
              if pink_recomb then pink_lallele{i} = blue_lallele{2};
              else                pink_lallele{i} = pink_lallele{2};
            end;
          end;
                                                                                                                         
          *-----------------------------------------------------------------------------------------;
          * Disease genotype index: 1=1/1, 2=1/2, 3=2/2, to feed into MUD{} array of genotypic means;
          *-----------------------------------------------------------------------------------------;
          dgene{i} = blue_dallele{i} + pink_dallele{i} - 1;
                                                                                                                         
          *-------------------------------------------------------------------;
          * Generate Phenotype from Disease Genotype for each person in family;
          *-------------------------------------------------------------------;
          call rannor(seed,x);
          phenos{i} = mud{dgene{i}} + de * x;

        end;
                                                                                                                         
        *-------------------------;
        * Output record (1/family);
        *-------------------------;
        output;
      end;
                                                                                                                         
      *--------------------------------;
      * Output last SEED for later use ;
      *--------------------------------;
      call symput('seed',trim(left(seed)));
      retain seed;
      drop i j;
    run;

                                                                                                                         
    proc print data=families(obs=10); title3 'Families'; run;

                                                                                                                       
    *------------------------------------------;
    * Now convert from Families to Individuals ;
    *------------------------------------------;
    data individs;
      set families;
                                                                                                                         
      *--------------------------------------------------------;
      * Copy the arrays defined before (above, in families dsn);
      *--------------------------------------------------------;
      array blue_dallele{*} blue_daf blue_dam blue_da1-blue_da&mkids;
      array pink_dallele{*} pink_daf pink_dam pink_da1-pink_da&mkids;
                                                                                                                         
      array blue_lallele{*} blue_laf blue_lam blue_la1-blue_la&mkids;
      array pink_lallele{*} pink_laf pink_lam pink_la1-pink_la&mkids;
                                                                                                                         
      length cdgf cdgm cdg1-cdg&mkids $10;
      length clgf clgm clg1-clg&mkids $10;
                                                                                                                         
      array cdgeno{*} cdgf cdgm cdg1-cdg&mkids;
      array clgeno{*} clgf clgm clg1-clg&mkids;
                                                                                                                         
      array phenos{*} pf  pm  p1-p&mkids;
                                                                                                                         
      do i=1 to &mkids+2;
        *----------------------------------------------------------------------------;
        * Convert the two numeric alleles into the character "x/y" genotype notation ;
        * for each of the disease gene and the linkage marker                        ;
        *----------------------------------------------------------------------------;
        cdgeno{i} = trim(left( min(of blue_dallele{i} pink_dallele{i}) ))||'/'||trim(left( max(of blue_dallele{i} pink_dallele{i}) ));
        clgeno{i} = trim(left( min(of blue_lallele{i} pink_lallele{i}) ))||'/'||trim(left( max(of blue_lallele{i} pink_lallele{i}) ));
        dgeno=cdgeno{i};
        lgeno=clgeno{i};
        QTp=phenos{i};
        QLp=(dgeno='2/2');


        output;
      end;
     run;
                                                                                                                         
     proc print data=individs; title3 'individs';
     run;
     proc export data=individs
     outfile="/home/ramua/project/linksim/linkdata.csv";
     run;
  

    data alleles;
      set families;
                                                                                                                         
      *--------------------------------------------------------;
      * Copy the arrays defined before (above, in families dsn);
      *--------------------------------------------------------;
      array blue_dallele{*} blue_daf blue_dam blue_da1-blue_da&mkids;
      array pink_dallele{*} pink_daf pink_dam pink_da1-pink_da&mkids;
                                                                                                                         
      array blue_lallele{*} blue_laf blue_lam blue_la1-blue_la&mkids;
      array pink_lallele{*} pink_laf pink_lam pink_la1-pink_la&mkids;
                                                                                                                         
      do i=1 to &mkids+2;
        *----------------------------------------------------------------------;
        * Output the two alleles (keeping the entire haplotype in the same obs);
        *----------------------------------------------------------------------;
        dallele = blue_dallele{i};
        lallele = blue_lallele{i};
        output;

        dallele = pink_dallele{i};
        lallele = pink_lallele{i};
        output;
      end;
     run;

     proc sort data=alleles; by i;
;
     OA
     *------------------------------;
     * Run simple Stat Analyses     ;
     *------------------------------;

     title3 'Is Disease Genotype Associated with Phenotype?';

     title4 'Histogram of QT PHENO';
     proc univariate data=individs plot;
      var QTp;
     run;

     proc univariate data=individs plot;
      class dgeno;
      var QTp;
     run;

     proc sort data=individs; by i pedid;


     title4 'GLM on Entire Dataset (ignoring family memberships)';
     proc glm data=individs;
       class dgeno;
       model QTp=dgeno;
     run;

     title4 'GLM on each member of Family';
     proc glm data=individs; by i;
       class dgeno;
       model QTp=dgeno;
     run;

     title4 'Mixed Model on Entire Dataset, accounting for Family memberships';
     proc mixed data=individs;
       class dgeno pedid;
       model QTp=dgeno;
       repeated  /subject=pedid type=un;
     run;


                                                                                  

     title3 'Is Linkage Genotype Associated with Phenotype?';
     title4 'GLM on Entire Dataset (ignoring family memberships)';
     proc glm data=individs;
       class lgeno;
       model QTp=lgeno;
     run;

     title4 'GLM on each member of Family';
     proc glm data=individs; by i;
       class lgeno;
       model QTp=lgeno;
     run;
     

     title4 'Mixed Model on Entire Dataset, accounting for Family memberships';
     proc mixed data=individs;
       class lgeno pedid;
       model QTp=lgeno;
       repeated  /subject=pedid type=un;
     run;




     title3 'Are the Disease Gene and Linkage Marker associated?';

     title4 'Contingency Table on Genotypes (ignoring family relationships)';
     proc freq data=individs;
       tables dgeno*lgeno/chisq;
     run;

     title4 'Contingency Table on Alleles (ignoring family memberships)';
     proc freq data=alleles;
       tables dallele*lallele/chisq;
     run;

     title4 'Contingency Table on each member of Family';
     proc freq data=alleles; by i;
       tables dallele*lallele/chisq;
     run;

     data alleles; set alleles;
       Y = dallele - 1;  * Y=0 or 1 indicating disease allele;
     run;

     proc GLIMMIX data = alleles;
       class pedid lallele;
       model Y(event=last) = lallele / dist=BINARY S OR;
       random _residual_ / subject=pedid type=un; 
       nloptions tech=none;
     run;

%mend link;

%link;
