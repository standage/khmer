//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include <math.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <string>

#include "khmer.hh"
#include "kmer_hash.hh"

using namespace std;

//
// _hash: hash a k-length DNA sequence into a 64-bit number.
//

namespace khmer
{

unsigned long long hash(const char * kmer, const WordLength k,
                   __uint128_t &_h, __uint128_t &_r)
{
    // sizeof(unsigned long long) * 8 bits / 2 bits/base
    if (!(k <= sizeof(unsigned long long)*4) || !(strlen(kmer) >= k)) {
      cout << "Testing Different K-Size" << endl; 
    }

    __uint128_t h = 0, r = 0;

    h |= twobit_repr(kmer[0]);
    r |= twobit_comp(kmer[k-1]);

    for (WordLength i = 1, j = k - 2; i < k; i++, j--) {
        h = h << 2;
        r = r << 2;

        h |= twobit_repr(kmer[i]);
        r |= twobit_comp(kmer[j]);
    }

    _h = h;
    _r = r;

   // hash both the forward and reverse bit representations
    hash<unsigned long long> kmer_ULL_hash;
    unsigned long long hashed_h = kmer_ULL_hash(h);
    unsigned long long hashed_r = kmer_ULL_hash(r);
    
    // return the lower hash value
    return uniqify_rc(hashed_h, hashed_r);
}

// _hash: return the maximum of the forward and reverse hash.

unsigned long long hash(const char * kmer, const WordLength k)
{
    unsigned long long h = 0;
    unsigned long long r = 0;

    return hash(kmer, k, h, r);
}

// _hash_forward: return the hash from the forward direction only.

unsigned long long hash_forward(const char * kmer, WordLength k)
{
    unsigned long long h = 0;
    unsigned long long r = 0;

    hash(kmer, k, h, r);
    return h;			// return forward only
}

//
// _revhash: given an unsigned int, return the associated k-mer.
//

std::string _revhash(unsigned long long hash, WordLength k)
{
    std::string s = "";
    return s;
}


};
