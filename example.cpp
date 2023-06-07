#include <iostream>
#include <vector>
#include "fft.cpp"

void print(std::vector<long long> a)
{
    for (long i = 0; i < a.size(); ++i)
        std::cout << a[i] << " ";
    std::cout << "\n";
}
int main()
{
    std::vector<long long> a{1, 2, 3}, b{1, 2, 3}; //for 1+2*x+3*x^2
    print(mult(a, b));//1 4 10 12 9 for 1+4*x+10*x^2+12*x^3+9*x^4
}
