class NucleotideHash()
{
  public:
    NucleotideHash(){
      for (int k = 0; k < hashvalues.length(); ++k) 
      {
        hashvalues[k] = value;
        value = dist(DRE);
      }
    }
    long* GetHashValue(int c) { return hashvalues[c]; }
    private:
      long mHashvalues = new [4];
}
