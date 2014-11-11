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

#include "khmer.hh"
#include "kmer_hash.hh"
#include "NucleotideHash.hh"

using namespace std;

//
// _hash: hash a k-length DNA sequence into a 64-bit number.
//

namespace khmer
{


HashIntoType _hash(const char * kmer, const WordLength k,
                   HashIntoType& _h, HashIntoType& _r)
{
  // constants FTW
  NucleotideHash *hasher = NucleotideHash();
  int n;
  // Using random for hash vals
  std::default_random_engine DRE;
  std::uniform_int_distribution<long> dist(1047483646,2147483646);

  // Instatiating the Nucleotide hash

  long getHashvalues(int c)
  {
    if(c > 3 || c < 0)
      throw khmer_exception();    
    return hasher.GetHashvalue(c);
  }// of getHashvalues

  void CyclicHash(int my_int)
  {
    n = my_int;
    if (n > WordLength) 
    {
      throw khmer_exception();
    }
  }

  // Bitwise Shifting 
  long fastLeftShiftN(long x)
  {
    return (x << n) | (x >> (WordLength - n));
  }

  long fastLeftShift1(long x)
  {
    return (x << 1) | (x >> (WordLength - 1));
  }

  long fastRightShiftN(long x)
  {
    return (x >> n) | (x << (WordLength - n));
  }

  long fastRightShift1(long x)
  {
    return (x >> 1) | (x << (WordLength - 1));
  }

  long eatRight(long hashval, int c)
  {
    hashval = fastLeftShift1(hashval);
    hashval ^= hasher.hashvalues[c];
    return hashval;
  }

  long eatLeft(long hashval, int c)
  {
    hashval = fastRightShift1(hashval ^ fastLeftShiftN(hasher.hashvalues[c]))
      return hashval;
  }

  long getInitialHashvalue(std::string s)
  {
    long hashvalue = 0;
    for (int i = 0; i < n; i++) 
    {
      hashvalue = fastLeftShift1(hashvalue);
      int idx = 0;
      switch (s.at(i)) 
      {
        case 'A':
        case 'a':
          idx = 0;
          break; 
        case 'C':
        case 'c':
          idx = 1;
          break; 
        case 'G':
        case 'g':
          idx = 2;
          break; 
        case 'T':
        case 't':
        case 'U':
        case 'u':
          idx = 3;
          break; 
        default:
          throw khmer_exception(); 
      }// of Switch
      hashvalue ^= hasher.getHashvalue(idx);
    }
    return hashvalue;
  }


  // sizeof(HashIntoType) * 8 bits / 2 bits/base
  if (!(k <= sizeof(HashIntoType)*4) || !(strlen(kmer) >= k)) {
    throw khmer_exception("Supplied kmer string doesn't match the underlying k-size.");
  }

  HashIntoType h = 0, r = 0;

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

  return uniqify_rc(h, r);
}

// _hash: return the maximum of the forward and reverse hash.

HashIntoType _hash(const char * kmer, const WordLength k)
{
  HashIntoType h = 0;
  HashIntoType r = 0;

  return khmer::_hash(kmer, k, h, r);
}

// _hash_forward: return the hash from the forward direction only.

HashIntoType _hash_forward(const char * kmer, WordLength k)
{
  HashIntoType h = 0;
  HashIntoType r = 0;


  khmer::_hash(kmer, k, h, r);
  return h;			// return forward only
}

//
// _revhash: given an unsigned int, return the associated k-mer.
//

std::string _revhash(HashIntoType hash, WordLength k)
{
  std::string s = "";

  unsigned int val = hash & 3;
  s += revtwobit_repr(val);

  for (WordLength i = 1; i < k; i++) {
    hash = hash >> 2;
    val = hash & 3;
    s += revtwobit_repr(val);
  }

  reverse(s.begin(), s.end());

  return s;
}


};
