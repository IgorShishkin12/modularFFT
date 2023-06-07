#include<vector>
long long mulmod(long long a, long long b, long long mod) {
    return (long long)((__int128(a) * __int128(b)) % __int128(mod));
}
std::vector<long long> mfft(std::vector<long long> dat, long long one, long long mod) {  // modular fast Fourier transform
    std::vector<long long> ans(dat.size());
    if (dat.size() <= 1) {//counting "by hand"
        ans[0] = 1;
        for (long i = 1; i < ans.size(); ++i)
            ans[i] = mulmod(ans[i - 1], one, mod);
        for (auto& i : ans) {
            long long ansval = 0;
            long long k = 1;
            for (auto j : dat) {
                ansval = (ansval + mulmod(k, j, mod)) % mod;
                k = mulmod(k, i, mod);
            }
            i = ansval;
        }
    } else {//fft step
        std::vector<long long> d1(dat.size() / 2), d2(dat.size() / 2);
        for (long i = 0; i < dat.size(); i += 2) {
            d1[i / 2] = dat[i];
            d2[i / 2] = dat[i + 1];
        }
        d1 = mfft(d1, mulmod(one, one, mod), mod);
        d2 = mfft(d2, mulmod(one, one, mod), mod);
        long long k = 1;
        for (long i = 0; i < ans.size(); ++i) {
            ans[i] = (d1[i % d1.size()] + mulmod(k, d2[i % d2.size()], mod)) % mod;
            k = mulmod(k, one, mod);
        }
    }
    return ans;
}
std::vector<long long> mult(std::vector<long long> a, std::vector<long long> b) {
    long long root = 9, invroot = 3502799710177053, module = 7881299347898369, degree = 49, degreeinv = 7881299347898355;  // settings

    long long len = max(a.size(), b.size());
    while (len & (len - 1))  // upper 2-th power
        len &= len - 1;
    long long lenpw = 2;
    while (len) {
        len /= 2;
        lenpw += 1;
    }
    //reducing degree to optimal
    while (degree > lenpw) {
        root = mulmod(root, root, module);
        invroot = mulmod(invroot, invroot, module);
        degree -= 1;
        degreeinv = mulmod(degreeinv, 2, module);
    }
    degree = 2LL << (degree - 1);
    a.resize(degree,0);
    b.resize(degree,0);
    //multiplication
    a = mfft(a, root, module);
    b = mfft(b, root, module);
    for (long i = 0; i < degree; ++i) {
        a[i] = mulmod(a[i], b[i], module);
    }
    a = mfft(a, invroot, module);
    for (auto& i : a) {
        i = mulmod(i, degreeinv, module);
    }
    //removing trailing zeros
    long long zeroCount = a.size()-1;
    while(zeroCount>=0&&a[zeroCount]==0){
        zeroCount--;
    }
    a.resize(zeroCount+1);
    if(a.empty())
        a.push_back(0);
    return a;
}
