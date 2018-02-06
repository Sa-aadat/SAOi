#############################################################################################
#                                                                                           #
#  Initialize                                                                               #
#                                                                                           #
#############################################################################################

rm -f *.out
make # clean

#############################################################################################
#                                                                                           #
#  The unconstrained test problems                                                          #
#                                                                                           #
#############################################################################################

# cp Examples/unconstrained/rosen1/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/quadratic/*.f      .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/quartic/*.f        .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/helical/*.f        .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/nlf3/*.f           .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/freudenroth/*.f    .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/freudenrothm/*.f   .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/powell_bad/*.f     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/brown_bad/*.f      .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/beale/*.f          .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unconstrained/wood/*.f           .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/rosenlarge/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi

#############################################################################################
#                                                                                           #
#  The unimodal test problems                                                               #
#                                                                                           #
#############################################################################################

cp Examples/unimodal/cam1/*.f                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/hock43/*.f              .
                                              rm -f Initialize.o Functions.o Gradients.o SAOi
                                              make
                                             ./SAOi
cp Examples/unimodal/cant1/*.f               .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/cant1a/*.f              .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/unimodal/cant2/*.f               .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unimodal/ccsa1/*.f               .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unimodal/ccsa2/*.f               .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
cp Examples/unimodal/default/*.f             .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/default_alt/*.f         .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/poly12/*.f              .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/rosen1/*.f              .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/unimodal/semiinf1/*.f            .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/unimodal/semiinf2/*.f            .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
cp Examples/unimodal/snake/*.f               .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/vdPlaats1/*.f           .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/vdPlaats2/*.f           .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/unimodal/vdPlaats3/*.f           .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/unimodal/weight1/*.f             .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi

#############################################################################################
#                                                                                           #
#  The multimodal test problems                                                             #
#                                                                                           #
#############################################################################################
# 
# cp Examples/multimodal/branin/*.f            .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/camel6/*.f            .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/goldprice/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/griewank2/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/griewank10/*.f        .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/hart3/*.f             .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/hart6/*.f             .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/partstamp1/*.f        .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/rastrigin/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/multimodal/schubert/*.f          .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
#############################################################################################
#                                                                                           #
#  The very large scale test problems                                                       #
#                                                                                           #
#############################################################################################
# 
# cp Examples/veryLarge/cam1/*.f               .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/cant2/*.f              .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/vdPlaats1/*.f          .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/vdPlaats2/*.f          .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/vdPlaats3/*.f          .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/weight1/*.f            .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/rosenlarge/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/veryLarge/springMass/*.f         .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
#############################################################################################
#                                                                                           #
#  The specialized test problems                                                            #
#                                                                                           #
#############################################################################################

cp Examples/specialized/dfree/bobyqa01/*.f   .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/dfree/newuoa01/*.f   .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/specialized/discont/sumsqr/*     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/specialized/discont/vdPlaats1/*  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/specialized/play/denseH/*        .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
cp Examples/specialized/play/qp1/*           .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/play/sala1/*         .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/specialized/nonconvex/ccsa1/*    .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# cp Examples/specialized/nonconvex/ccsa2/*    .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
cp Examples/specialized/sand/3bar/*          .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/sand/3barSAND/*      .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/sand/3barSANDsing/*  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/sand/3barSANDsing_vM/*  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

#############################################################################################
#                                                                                           #
#  The specialized separable test problems; exact Hessian info is given in diaHessUser.f    #
#                                                                                           #
#############################################################################################

cp Examples/specialized/separable/hock43/*   .
                                               rm -f Initialize.o Functions.o Gradients.o diaHessUser.o SAOi
                                               make
                                             ./SAOi
cp Examples/specialized/separable/cant1/*    .
                                               rm -f Initialize.o Functions.o Gradients.o diaHessUser.o SAOi
                                               make
                                             ./SAOi
rm diaHessUser.o
cp Repo/diaHess/diaHessUser.f                .

#############################################################################################
#                                                                                           #
#  The equality test problems                                                               #
#                                                                                           #
#############################################################################################

cp Examples/equalities/HS6/*                 .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS7/*                 .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS32/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS39/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS41/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS47/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS48/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS71/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS80/*                .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS107/*               .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS111/*               .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/HS112/*               .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/nonconvex1/*          .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/nonconvex2/*          .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
cp Examples/equalities/nonconvex3/*          .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
# cp Examples/equalities/vdplaats1Large/*      .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi

#############################################################################################
#                                                                                           #
#  The xAy test problems                                                                    #
#                                                                                           #
#############################################################################################
# 
# cp Examples/xAy/matstruc1/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# cp Examples/xAy/matstruc2/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# cp Examples/xAy/matstruc3/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# cp Examples/xAy/matstruc4/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/matstruc5/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/matstruc6/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# cp Examples/xAy/matstruc7/*                  .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# # cp Examples/xAy/testH0/*                     .
# #                                                rm -f Initialize.o Functions.o Gradients.o SAOi
# #                                                make
# #                                              ./SAOi

cp Examples/xAy/testH1/*                     .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

# cp Examples/xAy/testH2/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/testH3/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/testH4/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/testH5/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
#                                              
# cp Examples/xAy/testH6/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi
# 
# cp Examples/xAy/testH7/*                     .
#                                                rm -f Initialize.o Functions.o Gradients.o SAOi
#                                                make
#                                              ./SAOi

cat test.out
# diff test.out testref > diff.out
