This code implements(trying to) modularFFT(actually, polynomial multiplication, based on it), based on <a href = "https://www.csa.iisc.ac.in/~chandan/research/intMult.pdf">this article</a>.

You can use it as is, or with different parameters. Parameters you can find in file "Primes". Each line is a set of parameters, where Prime is prime of 1+j*2^i, root is root of one degree 2^i modulo Prime,1/root is inverted root modulo Prime. Set prime to module, root to root, 1/root to invroot, i to degree and 1/2^i to invdegree in settings in function mult.

Note: All input and output is modulo Prime!
