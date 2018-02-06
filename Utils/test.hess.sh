#############################################################################################
#                                                                                           #
#  Initialize                                                                               #
#                                                                                           #
#############################################################################################

rm -f *.out
make # clean


#############################################################################################
#                                                                                           #
#  The separable test problems - approximated                                               #
#                                                                                           #
#############################################################################################

cp Examples/unimodal/hock43/*.f               .
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
#rm diaHessUser.o
#cp Repo/diaHess/diaHessUser.f                .

cat test.out
# diff test.out testref > diff.out
