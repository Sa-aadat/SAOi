#############################################################################################
#                                                                                           #
#  Initialize                                                                               #
#                                                                                           #
#############################################################################################

rm -f *.out
make # clean

#############################################################################################
#                                                                                           #
#  The xAy test problems                                                                    #
#                                                                                           #
#############################################################################################

cp Examples/xAy/matstruc1/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/matstruc2/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/matstruc3/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/matstruc4/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
                                             
cp Examples/xAy/matstruc5/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi
                                             
cp Examples/xAy/matstruc6/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/matstruc7/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/testH0/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/testH1/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/testH2/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cp Examples/xAy/testH3/*                  .
                                               rm -f Initialize.o Functions.o Gradients.o SAOi
                                               make
                                             ./SAOi

cat test.out
# diff test.out testref > diff.out
